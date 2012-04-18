#! /usr/bin/env python
import multiprocessing

class ProcessConsumer(multiprocessing.Process):

    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                self.task_queue.task_done()
                break
            self.task_queue.task_done()
            # add the result of the Task's __call__
            # to the result queue
            self.result_queue.put(next_task())
        return


# class Task(object):
#     def __init__(self, var, func):
#         self.var = var
#     def __call__(self):
#         return '%s' % (self.var)
#     def __str__(self):
#         return '%s' % (self.var)
class Task(object):
    def __init__(self, func, args):
        self.func = func
        self.args = args
    def __call__(self):
        return self.func(*self.args)