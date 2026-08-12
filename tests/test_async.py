import asyncio
import subprocess

class Predictor:

    @staticmethod
    def _subprocess(cmd:list, shell:bool=False, tool_name:str=None) -> bool:
        try:
            if shell:
                cmd = ' '.join(cmd)
            res = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
            return res
        except Exception as e:
            print(f"ERROR: failed in {tool_name} execution, error={e}")
        return False
        
    @staticmethod
    async def f1(input, params):
        print("start f1")
        await asyncio.sleep(10)
        cmd = ['ls', ]
        Predictor._subprocess(cmd)
        return 'f1', input

    @staticmethod
    async def f2(input, *args):
        print("start f2")
        await asyncio.sleep(1)
        cmd = ['ls', ]
        Predictor._subprocess(cmd)
        return 'f2', input


async def run_worker(func, data_queue, result_queue, *args):
    while True:
        data = await data_queue.get()

        if data is None:
            data_queue.task_done()
            return
         
         # Send result to result queue
        try:
            result = await func(data, *args)
            await result_queue.put(result) 
        finally:
            data_queue.task_done()


async def process_results(result_queue, total_items, post_func=None):
    """
    Consumer the processes result in living time
    """
    processed = 0
    while processed < total_items:
        result = await result_queue.get()
        # Process immediately
        print(f'### {result}')
        if post_func:
            post_func(result)
        processed += 1
        
async def async_main(task_funcs:list, data:list, post_func=None, *args):
    result_queue = asyncio.Queue()
    
    # create input queue
    data_queue = asyncio.Queue()
    for item in data:
        await data_queue.put(item)

    # define worker tasks
    workers = []
    for func in task_funcs:
        task = asyncio.create_task(
            run_worker(func, data_queue, result_queue, *args)
        )
        workers.append(task)

    # Start result processor
    result_processor = asyncio.create_task(
        process_results(result_queue, len(data), post_func)
    )
    
    # Complete task: wait till all queued tasks done
    await data_queue.join()

    # complete and shutdown a worker if no data
    for _ in workers:
        await data_queue.put(None)

    # wait for all workers complete
    await asyncio.gather(*workers)

    # Wait for all results to be processed
    await result_processor

# wrap input data
data = ["A", "B", "C", "D", "E", "F", "G", "H"]
data = [{'id':i} for i in data]
# build funcs pool
funcs = [
    Predictor.f1,
    Predictor.f2,
]
params = {'host': 'am5'}
asyncio.run(async_main(funcs, data, None, params))
