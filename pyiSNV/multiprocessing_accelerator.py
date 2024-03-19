import multiprocessing
from math import ceil
from itertools import chain
import time
import os
from tqdm import tqdm

def factorial(n, **kwargs):
    # print(kwargs['extra_data'])

    if n < 0:
        raise ValueError("Factorial is not defined for negative numbers.")
    if n == 0:
        return 1
    result = 1
    for i in range(1, n + 1):
        result *= i

    return result

# def wrap_into_batch_func(batch, process_func):
#     # Print the current process ID
#     pid = f"Parent Process ID: {os.getppid()}, Current Process ID: {os.getpid()}"
#     print(pid)
#     res = []
#     # for i in tqdm(batch, desc=pid, miniters=10000):
#     #     res.append(process_func(i))
#     for i in batch:
#         res.append(process_func(i))
#     return res
#
#
# def chunkify(data, n):
#     """Yield n successive chunks from a list."""
#     size = ceil(len(data) / n)
#     for i in range(0, len(data), size):
#         yield data[i:i + size]
#
#
# def accelerate(data, process_func, num_processes):
#     """
#     Process the data in parallel while preserving order.
#     param data: An iterable of data to process.
#     return: A list of results in the same order as the input data.
#     """
#     print('begin multiprocessing')
#     start_time = time.time()
#
#     with multiprocessing.Pool(num_processes) as pool:
#         # Apply the function asynchronously and store the result objects
#         result_objects = [pool.apply_async(func=wrap_into_batch_func, args=(chunk, process_func))
#                           for chunk in chunkify(data, num_processes)]
#         # Retrieve the results in the order they were submitted
#         results = [result.get() for result in result_objects]
#         results = list(chain(*results))
#
#         # Wait for all tasks to complete
#         pool.close()
#         pool.join()
#
#     end_time = time.time()
#     runtime = end_time - start_time
#     print(f'{num_processes} processes: {runtime}')
#
#     return results



class AsyncMultiprocessingAccelerator():
    def __init__(self, ):
        # Don't encapsulate data into the class
        pass

    # ----------------------------------------------------------------
    # Overhead of Starting and Stopping Processes
    # Granularity of Tasks, Imbalance Workload, and Deadlocks
    # ----------------------------------------------------------------
    def accelerate(self, data, process_func, num_processes=20, **kwargs):

        print('begin multiprocessing')
        start_time = time.time()

        num_processes = num_processes or multiprocessing.cpu_count()
        with multiprocessing.Pool(num_processes) as pool:
            # Apply the function asynchronously and store the result objects
            result_objects = [pool.apply_async(func=self.wrap_into_batch_func, args=(chunk, process_func, i), kwds=kwargs)
                              for i, chunk in enumerate(self.chunkify_by_n(data, num_processes))]
            # Retrieve the results in the order they were submitted
            results = [result.get() for result in result_objects]
            results = list(chain(*results))

            # Wait for all tasks to complete
            pool.close()
            pool.join()

        end_time = time.time()
        runtime = end_time - start_time
        print(f'{num_processes} processes: {runtime}')

        return results

    @staticmethod
    def chunkify_by_n(data, n):
        """Yield n successive chunks from a list."""
        size = ceil(len(data) / n)
        for i in range(0, len(data), size):
            yield data[i:i + size]

    @ staticmethod
    def chunkify_by_size(data, size):
        """Yield n successive chunks from a list."""
        for i in range(0, len(data), size):
            yield data[i:i + size]

    @ staticmethod
    def wrap_into_batch_func(batch, process_func, id, **kwargs):
        # Print the current process ID
        pid = f"{id}. Parent Process ID: {os.getppid()}, Current Process ID: {os.getpid()} "
        res = []
        # print(pid)

        # for index, data in enumerate(batch):
        #
        #     res.append(process_func(data, **kwargs))
        #
        #     if index % 100 == 0:
        #         print(pid+f"{(index+1)/len(batch)*100:.0f}% {index+1}/{len(batch)}")

        for data in tqdm(batch, miniters=int(len(batch)*0.05) ,desc=pid):
            res.append(process_func(data, **kwargs))

        return res


if __name__ == "__main__":

    large_list = [200, 201, 202, 203, 204, 205, 206, 207, 208, 209] * 10_000_000

    #
    res = []
    start_time = time.time()
    for data in tqdm(large_list):
        res.append(factorial(data))
    end_time = time.time()
    runtime = end_time - start_time
    print(f'processes: {runtime}')
    print(res[0:10])


    # # ~25 minutes
    # res = []
    # start_time = time.time()
    # for index, data in enumerate(large_list):
    #
    #     res.append(factorial(data))
    #
    #     if index % 100000 == 0:
    #         print(f"{(index + 1) / len(large_list) * 100:.2f}% {index + 1}/{len(large_list)}")
    #
    # print(res[0:10])
    # end_time = time.time()
    # runtime = end_time - start_time
    # print(f'processes: {runtime}')


    # 116 seconds
    accelerator = AsyncMultiprocessingAccelerator()
    out = accelerator.accelerate(data=large_list, process_func=factorial, num_processes=40, extra_data='extra_data')
    print(out[0:10])

