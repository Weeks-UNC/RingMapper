

cdef struct READ:
    int start
    int stop
    int subcode
    char* read
    char* muts
    char* qual


cdef READ parseLine(char*, int)

cdef void fillReadMut(int[:], int[:], READ, int, int)

cdef void incrementArrays(int[:,::1], int[:,::1], int[:,::1], int[:], int[:])

