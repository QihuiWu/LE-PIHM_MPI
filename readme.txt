具体的代码修改如下：
1、在pihm.c文件的第78行插入“int rank, size;”。
2、在pihm.c文件的第79行插入“MPI_Status status;”。
3、在pihm.c文件的第80行插入“MPI_Init(&argc, &argv);”。
4、在pihm.c文件的第81行插入“MPI_Comm_rank(MPI_COMM_WORLD, &rank);”。
5、在pihm.c文件的第82行插入“MPI_Comm_size(MPI_COMM_WORLD, &size);”。
6、将pihm.c文件的第145行“sprintf(temp_path,"%s%d/",path,OutputID);”改成“sprintf(temp_path,"%s%d/",path,rank);”。
7、在pihm.c文件的第720行插入“MPI_Finalize();”。
8、在pihm.h文件的第32行插入“#include "mpi.h"”。
9、在read_alloc.c文件的第36行插入“int rank, size;”。
10、在read_alloc.c文件的第37行插入“MPI_Comm_rank(MPI_COMM_WORLD, &rank);”。
11、在read_alloc.c文件的第38行插入“MPI_Comm_size(MPI_COMM_WORLD, &size);”。
12、将read_alloc.c文件的第39行“int i, j;”改成“int i, j, k;”。
13、在read_alloc.c文件的第149行插入“for (k = 0; k <= rank; k++)”。
14、在read_alloc.c文件的第150行插入“{”。
15、在read_alloc.c文件的第166行插入“}”。
