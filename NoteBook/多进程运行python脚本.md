# temp

（1）使用multiprocessing进行多任务并行跑python脚本。  
```python
import multiprocessing

file_in = sys.argv[1]
all_cif_lst = get_lst(file_in)
random.shuffle(all_cif_lst)
with multiprocessing.Pool(processes=12) as pool:
    pool.map(pnprocess, all_cif_lst)  # 并行执行任务
```

（2）使用concurrent.futures进行多任务并行跑python脚本。感觉这个更好用一些。   
```python
import concurrent.futures
import time

# 定义一个任务函数
def task_function(n):
    print(f"Processing {n} on CPU...")
    time.sleep(2)  # 模拟耗时操作
    return n * n

# 主函数
if __name__ == "__main__":
    numbers = [1, 2, 3, 4, 5]  # 要处理的数据

    # 创建进程池执行器
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        # 使用 map 提交任务
        executor.map(task_function, numbers)
```