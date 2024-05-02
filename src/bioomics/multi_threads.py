import threading

def multi_threads(data_iter, func):
    threads = []
    for _data in data_iter:
        thread = threading.Thread(target=func, args=(_data,))
        threads.append(thread)
    
    for thread in threads:
        thread.start()
    
    for thread in threads:
        thread.join()

