# `bthread.h` - Usage Guide

> Lightweight worker-thread pool helpers implemented as a set of macros and small inline routines.
>
> This document explains what the header does, when to use it, and how to use it with examples and notes about portability and threading semantics.

---

## What is this library

The header provides a macro-based way to create simple worker thread pools with minimal boilerplate. It exposes a pattern that:

* Generates a typed worker thread structure and pool for a user-specified argument type.
* Uses C11 atomics for lightweight state signalling between the application and worker threads.
* Optionally binds worker threads to NUMA nodes when `<numa.h>` is available.

---

## When to use it

Use this header when you need:

* A small, fast, simple thread pool where each worker is persistent and waits for work signalled by an atomic variable.
* Low-latency coordination between a producer (main thread) that triggers work on specific workers and the workers themselves.
* Per-worker arguments that are known at pool creation time.
* Optional NUMA affinity (if available on the system) to bind workers to nodes.

Avoid using this header when:

* You need a dynamic work-queue where tasks are enqueued and workers pull from the queue.
* You require advanced features like work-stealing, dynamic scaling of worker count, cancellation-of-in-flight-work, or thread pool scheduling policies.
* You prefer blocking synchronization primitives (mutex+condvar) rather than spin-waiting — spin-waiting wastes CPU if the wait is long.

---

## Requirements

* C11 compiler (for `stdatomic.h`).
* `pthread` library.
* Optional: `libnuma` and `<numa.h>` for NUMA binding.

Compile/link flags common examples:

```sh
# Without NUMA:
cc -std=c11 -pthread -o myprog myprog.c

# If using libnuma (and header is present):
cc -std=c11 -pthread -lnuma -o myprog myprog.c
```

---

## Macro API overview

The header exposes the following macro and generated symbols:

### `INIT_BTHREAD(name, type)`

Generates the following types and functions (where `name` is your pool name, `type` is the pointer argument type used by worker functions):

* `void name##_setarg(name##_bpool_t *btp, size_t i, bthread_func_t func, type *arg);`
* `type *name##_getarg(name##_bpool_t *btp, size_t i);`
* `void name##_start(name##_bpool_t *btp, size_t i);` - signal single worker to run.
* `void name##_wait(name##_bpool_t *btp, size_t i);` - spin until worker reports DONE, then resets to IDLE.
* `name##_bpool_t *name##_init(bthread_func_t func, type **args, size_t thread_cnt);` - create pool and spawn threads.
* `void name##_destroy(name##_bpool_t *btp);` - stop workers and free resources.
---

## Important implementation notes and gotchas

* **Spin-waiting:** `name##_wait()` and the worker loop use busy-wait loops with `cpu_relax()` — this is low-latency but can be CPU intensive if workers spend long times in WAIT. Use only for short, frequent tasks or on dedicated cores.
* **Atomic memory ordering:** Status transitions use a mix of `memory_order_acquire`, `memory_order_release`, and `memory_order_acq_rel` for correctness of the simple state machine.
* **Minimum thread count:** `__MINIMUM_THREAD_COUNT` is `1`. The `init` function enforces at least one worker.
* **Worker args array:** `name##_init` expects an array of pointers `type **args` of length `thread_cnt`. If you store shared state, ensure it is valid for the lifetime of the workers.

---

## Example usage

This example demonstrates a pool that runs a per-thread function that prints a message and the thread id.

```c
#include "bi_thread.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

typedef struct {
    int id;
    const char *msg;
} myarg_t;

/* instantiate pool types and functions for myarg_t */
INIT_BTHREAD(my, myarg_t)

/* worker function signature required: void (*)(void *arg) */
void worker_fn(void *arg_void) {
    myarg_t *arg = (myarg_t*)arg_void;
    // Do some work
    printf("Worker %d: %s\n", arg->id, arg->msg);
    // simulate work
    sleep(1); // avoid long sleeps for spinwait pools on production
}

int main(void) {
    size_t n = 4;
    myarg_t *args[4];

    for (size_t i = 0; i < n; ++i) {
        args[i] = malloc(sizeof(myarg_t));
        args[i]->id = (int)i;
        args[i]->msg = "hello from worker";
    }

    // initialize pool, passing the same function pointer and per-worker args
    my_bpool_t *pool = my_init(worker_fn, (myarg_t**)args, n);
    if (!pool) {
        fprintf(stderr, "failed to init pool\n");
        return 1;
    }

    // start each worker and wait for completion
    for (size_t i = 0; i < n; ++i) {
        my_start(pool, i);   // signal worker i to run
    }

    for (size_t i = 0; i < n; ++i) {
        my_wait(pool, i);    // wait until worker i finished (and resets to IDLE)
    }

    // reuse workers: change args or call my_setarg() and start again as needed.

    my_destroy(pool);

    for (size_t i = 0; i < n; ++i) free(args[i]);

    return 0;
}
```

Compile:

```sh
cc -std=c11 -pthread -o myprog myprog.c
./myprog
```

---

## Summary

This header provides a practical, low-latency per-worker signalling thread pool implementation, ideal for small tasks or dedicated core setups where tasks are short and frequent. It trades CPU usage (spin-wait) for low wake-up latency and simplicity. Use it where simplicity and performance matter, but avoid it when you need general-purpose task queues or power-efficient waits.