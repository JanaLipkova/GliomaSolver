#ifndef _NULL_MPI_H_
#define _NULL_MPI_H_

#define NULL_MPI	1

/*** Public ***/
typedef int MPI_Datatype;
#define MPI_CHAR           ((MPI_Datatype)0x4c000101)
#define MPI_SIGNED_CHAR    ((MPI_Datatype)0x4c000118)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)0x4c000102)
#define MPI_BYTE           ((MPI_Datatype)0x4c00010d)
#define MPI_WCHAR          ((MPI_Datatype)0x4c00040e)
#define MPI_SHORT          ((MPI_Datatype)0x4c000203)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)0x4c000204)
#define MPI_INT            ((MPI_Datatype)0x4c000405)
#define MPI_UNSIGNED       ((MPI_Datatype)0x4c000406)
#define MPI_LONG           ((MPI_Datatype)0x4c000407)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)0x4c000408)
#define MPI_FLOAT          ((MPI_Datatype)0x4c00040a)
#define MPI_DOUBLE         ((MPI_Datatype)0x4c00080b)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)0x4c000c0c)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)0x4c000809)
#define MPI_UNSIGNED_LONG_LONG ((MPI_Datatype)0x4c000819)
#define MPI_LONG_LONG      MPI_LONG_LONG_INT

#define MPI_PACKED         ((MPI_Datatype)0x4c00010f)
#define MPI_LB             ((MPI_Datatype)0x4c000010)
#define MPI_UB             ((MPI_Datatype)0x4c000011)

/*
   The layouts for the types MPI_DOUBLE_INT etc are simply
   struct {
       double var;
       int    loc;
   }
   This is documented in the man pages on the various datatypes.
 */
#define MPI_FLOAT_INT         ((MPI_Datatype)0x8c000000)
#define MPI_DOUBLE_INT        ((MPI_Datatype)0x8c000001)
#define MPI_LONG_INT          ((MPI_Datatype)0x8c000002)
#define MPI_SHORT_INT         ((MPI_Datatype)0x8c000003)
#define MPI_2INT              ((MPI_Datatype)0x4c000816)
#define MPI_LONG_DOUBLE_INT   ((MPI_Datatype)0x8c000004)


/* Fortran types */
#define MPI_COMPLEX           ((MPI_Datatype)1275070494)
#define MPI_DOUBLE_COMPLEX    ((MPI_Datatype)1275072546)
#define MPI_LOGICAL           ((MPI_Datatype)1275069469)
#define MPI_REAL              ((MPI_Datatype)1275069468)
#define MPI_DOUBLE_PRECISION  ((MPI_Datatype)1275070495)
#define MPI_INTEGER           ((MPI_Datatype)1275069467)
#define MPI_2INTEGER          ((MPI_Datatype)1275070496)
#define MPI_2COMPLEX          ((MPI_Datatype)1275072548)
#define MPI_2DOUBLE_COMPLEX   ((MPI_Datatype)1275076645)
#define MPI_2REAL             ((MPI_Datatype)1275070497)
#define MPI_2DOUBLE_PRECISION ((MPI_Datatype)1275072547)
#define MPI_CHARACTER         ((MPI_Datatype)1275068698)

/* Size-specific types (see MPI-2, 10.2.5) */
#define MPI_REAL4             ((MPI_Datatype)0x4c000427)
#define MPI_REAL8             ((MPI_Datatype)0x4c000829)
#define MPI_REAL16            ((MPI_Datatype)MPI_DATATYPE_NULL)
#define MPI_COMPLEX8          ((MPI_Datatype)0x4c000828)
#define MPI_COMPLEX16         ((MPI_Datatype)0x4c00102a)
#define MPI_COMPLEX32         ((MPI_Datatype)MPI_DATATYPE_NULL)
#define MPI_INTEGER1          ((MPI_Datatype)0x4c00012d)
#define MPI_INTEGER2          ((MPI_Datatype)0x4c00022f)
#define MPI_INTEGER4          ((MPI_Datatype)0x4c000430)
#define MPI_INTEGER8          ((MPI_Datatype)0x4c000831)
#define MPI_INTEGER16         ((MPI_Datatype)MPI_DATATYPE_NULL)

/*** Private ***/


/* Communicators */
typedef int MPI_Comm;
#define MPI_COMM_NULL  ((MPI_Comm)0x04000000)
#define MPI_COMM_WORLD ((MPI_Comm)0x44000000)
#define MPI_COMM_SELF ((MPI_Comm)0x44000001)

/* for info */
typedef int MPI_Info;
#define MPI_INFO_NULL         ((MPI_Info)0x1c000000)

/* The order of these elements must match that in mpif.h */
typedef struct MPI_Status {
    int count;
    int cancelled;
    int MPI_SOURCE;
    int MPI_TAG;
    int MPI_ERROR;

} MPI_Status;


/* RMA and Windows */
typedef int MPI_Win;
#define MPI_WIN_NULL ((MPI_Win)0x20000000)


/* MPI's error classes */
#define MPI_SUCCESS          0      /* Successful return */

/* Definitions that are determined by configure. */          
typedef int MPI_Fint;

/* Handle conversion types/functions */   

/* Programs that need to convert types used in MPICH should use these */
#define MPI_Type_f2c(datatype) (MPI_Datatype)(datatype)

#define MPI_ARGV_NULL (char **)0       

/* For supported thread levels */
#define MPI_THREAD_SINGLE 0
#define MPI_THREAD_FUNNELED 1
#define MPI_THREAD_SERIALIZED 2
#define MPI_THREAD_MULTIPLE 3




#define MPI_ANY_SOURCE  (-2)
#define MPI_LOCK_EXCLUSIVE  234
#define MPI_LOCK_SHARED     235


#define MPID_Datatype_get_basic_size(a) (((a)&0x0000ff00)>>8)

#define MPI_Send(a1,a2,a3,a4,a5,a6)
#define MPI_Ssend(a1,a2,a3,a4,a5,a6)
#define MPI_Recv(a1,a2,a3,a4,a5,a6,a7)	((*a7).MPI_ERROR = MPI_SUCCESS)
#define MPI_Type_size(type, size)	do {*size = MPID_Datatype_get_basic_size(type);} while(0)
#define MPI_Barrier(a1)
#define MPI_Bcast(a1,a2,a3,a4,a5)
#define MPI_Comm_size(a1,a2)		do {*a2 = 1;} while (0)
#define MPI_Comm_rank(a1,a2)		do {*a2 = 0;} while (0)
#define MPI_Intercomm_merge(a1,a2,a3)	do {} while (0)

#define MPI_Get_processor_name(a1,a2)	gethostname(a1, a2)
#define MPI_Finalize()
#define MPI_Abort(a1,a2)

#define MPI_Finalized(a1)		do {*a1 = 0;} while (0)
#define MPI_Initialized(a1)		do {*a1 = 1;} while (0)
#define MPI_Query_thread(a1)		do {*a1 = MPI_THREAD_MULTIPLE;} while (0)
#define MPI_Comm_get_parent(a1)		do {*a1 = MPI_COMM_NULL; } while (0)
   

#define MPI_Get(a1,a2,a3,a4,a5,a6,a7,a8)
#define MPI_Win_create(a1,a2,a3,a4,a5,a6)
#define MPI_Win_lock(a1,a2,a3,a4)
#define MPI_Win_unlock(a1,a2)

#define MPI_Alloc_mem(a1,a2,a3)
#define MPI_Init_thread(a1,a2,a3,a4)	do {*a4 = a3;} while(0)
#define MPI_Comm_dup(a1,a2)



#if 1
typedef int MPI_Request;

#define MPI_Irecv(a1,a2,a3,a4,a5,a6,a7)		1
#define MPI_Test(a1,a2,a3)

#endif

#if 1
//typedef int MPI_Fint;
#define MPI_DATATYPE_NULL  ((MPI_Datatype)0x0c000000)
#define MPI_Init(a1,a2)	
#define MPI_MAX_PROCESSOR_NAME 128      
#define MPI_Reduce(a1,a2,a3,a4,a5,a6,a7)

typedef int MPI_Op;
#define MPI_SUM     (MPI_Op)(0x58000003)

double MPI_Wtime(void);

#define MPI_Allgather(a1,a2,a3,a4,a5,a6,a7)

#define MPI_UNIVERSE_SIZE    0x64400009

#endif
#endif
