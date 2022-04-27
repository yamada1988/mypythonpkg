import threading
from threading import Thread

# add stop(), returnvalue to threading.Thread class.
class MyThread(Thread):
    def __init__(self, group=None, target=None, name=None,
                 args=(), kwargs={}, Verbose=None):
        Thread.__init__(self, group, target, name, args, kwargs, Verbose)
        self._return = None
        self.stop_event = threading.Event()
        self.setDaemon(True)

    def run(self):
        if self._Thread__target is not None:
            self._return = self._Thread__target(*self._Thread__args,
                                                **self._Thread__kwargs)
    def join(self):
        Thread.join(self)
        return self._return

    def stop(self):
        self.stop_event.set()
