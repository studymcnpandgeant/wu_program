wu_program
==========  
![smile](smile.png)  
这里主要存储一些前期的PETSc小程序，关键性的测试算例等。  

# 2018.10.10  
上传了4个程序
1. 2群中子稳态计算程序，均匀物质
2. 热工部分耦合计算算例
3. slepc中子本征值计算程序
4. PETSc矩阵读取例程  
这里需要注意的是：
注意每个文件的输入命令
 `test args: `

# 2018.10.11
学习了2个TS算例，对TS Object有一个大体的了解。  
上传了一个小程序，可以支持Matrix-Free以及Field Split PC  
在我的理解中就是JFNK方法再加上块Jacobi预处理  
可以借鉴该算例，修改耦合计算程序


# 2018.10.12
多文件C++程序处理  
目前效果不是很好，仅能实现把部分程序放入一个头文件中。  
可以参考  
[Multifiles](http://www.itkeyword.com/doc/5923142964970566764/creating-a-library-file-in-makefile-and-compiling-after-that)  
上述文档中说明了，如何把以及编译好的库文件链接到PETSc文件中。  
现在的问题是如何把自己的类或者函数封装为库文件，尤其是含有某些PETSc对象的函数，如何生成.o并生成lib？


# 2018.10.16
中子部分并行程序的修改，主要体现在以下不同：  
* keff变成一个物理场，但是由于这个物理场残差函数是每个网格都相同的，所以可以保证是一个均匀的物理场。
* 把FormFuntion变成区域计算，在每个处理器上只计算属于该处理器上的内容，这有可能涉及到映射点的问题。
* 程序中使用了DMCompositeScatter来获取局部向量，但是官方建议是使用`DMGlobalToLocal()` and `DMLocalToGlobal()` 来获得映射点。
* 有一个致命的问题：在计算keff的时候需要每个网格的信息，这是仅靠映射点以及不够了，应该如何处理？



# 2018.10.17
对指针的理解：  
`2164: PetscErrorCode  VecGetArray2d(Vec x,PetscInt m,PetscInt n,PetscInt mstart,PetscInt nstart,PetscScalar **a[])`  
`2165: { `   
`2167:   PetscInt       i,N;`      
`2168:   PetscScalar    *aa;`   
`2174:   VecGetLocalSize(x,&N);`  
`2175:   if (m*n != N) SETERRQ3(PETSC_COMM_SELF,PETSC_ERR_ARG_INCOMP,"Local array size %D does not match 2d array dimensions %D by %D",N,m,n);  `  
`2176:   VecGetArray(x,&aa); `   
`2178:   PetscMalloc1(m,a);`    
`2179:   for (i=0; i<m; i++) (*a)[i] = aa + i*n - nstart; `   
`2180:   *a -= mstart;`    
`2181:   return(0); `   
`2182: }`  
这里的`**a[]`应该指向一个二维数组，从声明上可以看出这是一个数组，每个元素是指针的指针，从函数里可以看出这个数组维度是m。  
`2179`行表示的是，对每个指针的起始位置的地址赋予初值。  


# 2018.10.18  
1. 上传了两个小程序，用来读取SN2D输入卡中的数据，并赋值给结构体数组中。
2. 学习指针  
2.1  指针名表示地址，`*`运算符被称为间接值或者解除引用算符。  
2.2  `int* ptr`强调指向Int的一个值，该值的地址为ptr。  
2.3  使用`new`来分配内存  
  在运行阶段分配未命名的内存以存储值  
  `typeName* pointer_name = new typeName`  
  eg.`int* pt = new int`  
   `*pt = 1000`  
   先分配了内存地址，再使用该内存地址的指针指向一个数据。  
 2.4  指针数组运算符  
 `arrayname[i]  becomes *(arrayname + i)`  
 2.5  如果结构标识符是结构名，则使用句点运算符；如果标识符是指向结构的指针，则使用箭头标识符。  
3. 反思算例/snes/ex28.c  
程序中最初定义了一个指向结构体UserCtx的指针`typedef struct _UserCtx* User`  
这里定义了一个指向结构体的指针类型，这个类型之后在主函数中声明  
`User user`这里的user是一个指针（地址），这个指针指向结构体_UserCtx  
指针定义了，但是需要对这个指针初始化，使它指向一个实际的结构体 **allocate space for struct**  
`PetscNew(&user)`作用类似于new，给user这个地址分配一个内存空间。  
由于user这里是地址（指针），所以使用结构体内变量使用箭头标识符。  


# 2018.10.19
对引用的理解：  
