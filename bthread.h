#ifndef __BI_THREAD__
#define __BI_THREAD__

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>
#include <pthread.h>
#include <stdatomic.h>

#define __WAIT_TIME_LITTLE      256
#define __WAIT_TIME_SHORT       4096 
#define __WAIT_TIME_MODERATE    32768
#define __MINIMUM_THREAD_COUNT  1

#if __has_include(<numa.h>)
#  include <numa.h>
#  define HAVE_NUMA_H 1
#else
#  define HAVE_NUMA_H 0
#endif

static inline void cpu_relax(void) {
    for (int i = 0; i < __WAIT_TIME_LITTLE; i++)
#if defined(__x86_64__)
        __asm__ __volatile__("pause" ::: "memory");
#elif defined(__aarch64__)
        __asm__ __volatile__("yield" ::: "memory");
#else
        __asm__ __volatile__("" ::: "memory");
#endif
}

#if HAVE_NUMA_H
#define bind_to_node(node_id) {                                                                                                 \
    struct bitmask *bm = numa_allocate_nodemask();                                                                              \
    numa_bitmask_setbit(bm, node_id % (numa_max_node()+1));                                                                     \
    numa_bind(bm);                                                                                                              \
    numa_free_nodemask(bm);                                                                                                     \
}
#else
#define bind_to_node(node_id) ((void)(node_id))
#endif


typedef void (*bthread_func_t)(void *arg);

#define INIT_BTHREAD(name, type)                                                                                                \
                                                                                                                                \
    enum {                                                                                                                      \
        __BTHREAD_IDLE   = 0,                                                                                                   \
        __BTHREAD_RUN    = 1,                                                                                                   \
        __BTHREAD_DONE   = 2,                                                                                                   \
        __BTHREAD_WORK   = 3,                                                                                                   \
        __BTHREAD_WAIT   = 4,                                                                                                   \
        __BTHREAD_EXIT   = 5,                                                                                                   \
    };                                                                                                                          \
                                                                                                                                \
    typedef struct {                                                                                                            \
        bthread_func_t  func;                                                                                                   \
        type *          arg;                                                                                                    \
        size_t          thd_id;                                                                                                 \
        atomic_int      status;                                                                                                 \
        pthread_t       thread;                                                                                                 \
    } name##_bthread_t;                                                                                                         \
                                                                                                                                \
    typedef struct {                                                                                                            \
        name##_bthread_t *workers;                                                                                              \
        size_t thread_cnt;                                                                                                      \
    } name##_bpool_t;                                                                                                           \
                                                                                                                                \
    static inline void name##_setarg(name##_bpool_t *btp, size_t i, bthread_func_t func, type *arg) {                           \
        if (!btp || i >= btp->thread_cnt) return;                                                                               \
        btp->workers[i].func = func;                                                                                            \
        btp->workers[i].arg = arg;                                                                                              \
        btp->workers[i].thd_id = i;                                                                                             \
        atomic_init(&(btp->workers[i].status), __BTHREAD_IDLE);                                                                 \
    }                                                                                                                           \
                                                                                                                                \
    static inline type *name##_getarg(name##_bpool_t *btp, size_t i) {                                                          \
        if (!btp || i >= btp->thread_cnt) return NULL;                                                                          \
        return btp->workers[i].arg;                                                                                             \
    }                                                                                                                           \
                                                                                                                                \
    static inline void name##_start(name##_bpool_t *btp, size_t i) {                                                            \
        if (!btp || i >= btp->thread_cnt) return;                                                                               \
        atomic_store_explicit(&(btp->workers[i].status), __BTHREAD_RUN, memory_order_release);                                  \
    }                                                                                                                           \
                                                                                                                                \
    static inline void name##_wait(name##_bpool_t *btp, size_t i) {                                                             \
        if (!btp || i >= btp->thread_cnt) return;                                                                               \
        name##_bthread_t *worker = &btp->workers[i];                                                                            \
        while (atomic_load_explicit(&worker->status, memory_order_acquire) != __BTHREAD_DONE)                                   \
            cpu_relax();                                                                                                        \
        /* reset to idle */                                                                                                     \
        atomic_store_explicit(&worker->status, __BTHREAD_IDLE, memory_order_release);                                           \
    }                                                                                                                           \
                                                                                                                                \
    void *name##__worker(void *arg) {                                                                                           \
        name##_bthread_t *worker = (name##_bthread_t *)arg;                                                                     \
        bind_to_node(worker->thd_id);                                                                                           \
                                                                                                                                \
        for ( ;; ) {                                                                                                            \
            int s;                                                                                                              \
            while ((s = atomic_load_explicit(&worker->status, memory_order_acquire)) != __BTHREAD_RUN && s != __BTHREAD_EXIT)   \
                cpu_relax();                                                                                                    \
                                                                                                                                \
            if (s == __BTHREAD_EXIT) break;                                                                                     \
                                                                                                                                \
            int expected = __BTHREAD_RUN;                                                                                       \
            if (!atomic_compare_exchange_strong_explicit(                                                                       \
                    &worker->status,     /* pointer to atomic variable */                                                       \
                    &expected,           /* expected value (input/output) */                                                    \
                    __BTHREAD_WORK,      /* value to store if match */                                                          \
                    memory_order_acq_rel,/* success memory order */                                                             \
                    memory_order_acquire /* failure memory order */                                                             \
            )) continue;                                                                                                        \
                                                                                                                                \
            /* do the job */                                                                                                    \
            worker->func(worker->arg);                                                                                          \
                                                                                                                                \
            /* mark finished */                                                                                                 \
            atomic_store_explicit(&worker->status, __BTHREAD_DONE, memory_order_release);                                       \
        }                                                                                                                       \
        return NULL;                                                                                                            \
    }                                                                                                                           \
                                                                                                                                \
    static inline name##_bpool_t *name##_init(bthread_func_t func, type **args, size_t thread_cnt) {                            \
        if (!func || thread_cnt <= 0) return NULL;                                                                              \
        name##_bpool_t *btp = calloc(1, sizeof(*btp));                                                                          \
        if (!btp) return NULL;                                                                                                  \
        thread_cnt = (thread_cnt < __MINIMUM_THREAD_COUNT) ? __MINIMUM_THREAD_COUNT : thread_cnt;                               \
        btp->thread_cnt = thread_cnt;                                                                                           \
        btp->workers = calloc(thread_cnt, sizeof(name##_bthread_t));                                                            \
        if (!btp->workers) { free(btp); return NULL; }                                                                          \
                                                                                                                                \
        for (size_t i = 0; i < thread_cnt; i++) {                                                                               \
            name##_setarg(btp, i, func, args[i]);                                                                               \
            pthread_create(&(btp->workers[i].thread), NULL, name##__worker, &(btp->workers[i]));                                \
        }                                                                                                                       \
        return btp;                                                                                                             \
    }                                                                                                                           \
                                                                                                                                \
    static inline void name##_destroy(name##_bpool_t *btp) {                                                                    \
        if (!btp) return;                                                                                                       \
                                                                                                                                \
        for (int i = 0; i < btp->thread_cnt; i++) {                                                                             \
            name##_bthread_t *worker = &btp->workers[i];                                                                        \
            atomic_store_explicit(&worker->status, __BTHREAD_EXIT, memory_order_release);                                       \
        }                                                                                                                       \
        for (size_t i = 0; i < btp->thread_cnt; i++)                                                                            \
            pthread_join(btp->workers[i].thread, NULL);                                                                         \
                                                                                                                                \
        free(btp->workers);                                                                                                     \
        free(btp);                                                                                                              \
    }


#ifdef __cplusplus
}
#endif

#endif