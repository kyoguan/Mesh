Mesh: Compacting Memory Management for C/C++
============================================

原版的Mesh可以访问[https://github.com/plasma-umass/Mesh]

Latest release tag: v1.0.3 (2021-02-26).

Release
------------------
* 2021-02-26 `v1.0.3` mesh一次只对一个class, 减少了mesh消耗的cpu.
* 2021-02-18 `v1.0.0` 新增实现了fork之后的COW。

这里做了几点对比原作的改进

1. 修正了很多bug使得他可以通过各个内存管理器的压测程序，其中严重的bug（会导致crash）已经Merge回原版。
2. 大幅的性能改进，在实际应用的时候性能已经跟jemalloc相当，某些程序中还要好于jemalloc
3. 增加了一个可以用线程辅助内存整理的选项，可以通过环境变量打开，线程模式可以令性能大幅提升，但会额外多消耗1%左右的cpu

```
MESH_FREEPHYS_THREAD=1  LD_PRELOAD=~/work/Mesh/build/lib/libmesh.so ./test-stress 64 10 20000000
```

4. 在MESH_BACKGROUND_THREAD选项下，现在可以输出更多的内存信息，通过向程序发送SIGUSR1 和 SIGUSR2 得到不同等级的信息

```
MESH_BACKGROUND_THREAD=1 MESH_FREEPHYS_THREAD=1 MESH_PERIOD_MS=100 LD_PRELOAD=~/work/Mesh/build/lib/libmesh.so ./test-stress 64 10 20000000


MESH_BACKGROUND_THREAD=1 MESH_FREEPHYS_THREAD=1 LD_PRELOAD=~/work/Mesh/build/lib/libmesh.so ./test-stress 64 10 20000000
Using 64 threads with a 10% load-per-thread and 20000000 iterations
- iterations left: 19999990,  use - 0.122959
- iterations left: 19999980,  use - 0.138058
- iterations left: 19999970,  use - 0.136929
- iterations left: 19999960,  use - 0.114705
- iterations left: 19999770,  use - 0.125529
MeshInfo ++++++++++ partial class:1 , length:36
MeshInfo ++++++++++ partial class:2 , length:24
MeshInfo ++++++++++ partial class:4 , length:36
MeshInfo ++++++++++ partial class:8 , length:67
MeshInfo ++++++++++ partial class:19, length:1
MeshInfo ++++++++++ partial class:23, length:1
MeshInfo ++++++++++ partial class:27, length:1
MeshInfo ++++++++++ partial class:31, length:3
MeshInfo ++++++++++ partial class:35, length:4
```

5. 内存管理块的大小有变动，使得内存块能100%被利用，原版由于块的大小在某些时候不整除，导致有部分浪费，这部分的变动来自于参考jemalloc的设置

6. 改进span的合并和扫描方式，现在可以分步扫描，更加平滑，更少的lock
7. 改进了VM地址的回收，Mesh后的虚地址可以更快的提前回收，可以减少整体vm地址的增长
8. 改进了fork的速度，大概可以有8倍速度的提升


How to build?
------------------


```
git clone https://github.com/kyoguan/Mesh
cd Mesh
cmake .
make -j
```


注意
------------------

1. Mesh不能跟THP一起用，所以如果你系统打开了THP，请关掉，或者设置为madvise模式
2. 由于使用Mesh的程序需要消耗大量maps的结构，所以需要调大 vm.max_map_count
```
sudo sysctl -w vm.max_map_count=100000000
```


为什么
------------------

1. 为什么Mesh一旦跑起来，VIRT就很大

```
%Cpu(s): 25.8 us, 25.8 sy,  0.0 ni, 48.4 id,  0.0 wa,  0.0 hi,  0.0 si,  0.0 st
MiB Mem :  62407.7 total,  56860.3 free,   1144.0 used,   4403.5 buff/cache
MiB Swap:   2048.0 total,   2048.0 free,      0.0 used.  60394.3 avail Mem

 进程号 USER      PR  NI    VIRT    RES    SHR    %CPU  %MEM     TIME+ COMMAND
  65421 kyo       20   0  134.8g 175064 168724 S 400.0   0.3   0:14.25 test-stress
```

  这个不用担心，因为这个只是虚地址的空间，它实际上不占用真实内存


2. 为什么对比其他分配器，RES并没有明显区别？

  这里有个知识点，RSS跟PSS的区别。RSS是代表有多少虚拟地址实际占用了物理内存，PSS是表示实际占用了多少物理内存。举个例子：如果32K虚拟地址映射到 16K的物理内存上，那么RSS就是32K, PSS就是16K。遗憾的是系统没有很方便统计PSS的工具，只有统计RSS。

  ```
MALLOCSTATS=2 LD_PRELOAD=~/work/Mesh/build/lib/libmesh.so ./test-stress
Using 32 threads with a 10% load-per-thread and 50 iterations
- iterations left:  40,  use - 0.084576
- iterations left:  30,  use - 0.090809
- iterations left:  20,  use - 0.080349
- iterations left:  10,  use - 0.082333
- iterations left:   0,  use - 0.093384
MESH COUNT:         3920
Meshed MB (total):  15.3
Meshed pages HWM:   243
Meshed MB HWM:      0.9
MH Alloc Count:     117810
MH Free  Count:     498610
MH High Water Mark: 10311
  ```

  可以通过MALLOCSTATS环境变量，让程序输出统计信息，这里显示的就是15MB是合并过的，所以实际是不占物理内存的。经验数据上来说，Mesh实际可以减少大概10%~30%的物理内存占用，尤其对于长期运行的程序，因为他的整理过程是一直持续的，达到不断整理碎片的作用。