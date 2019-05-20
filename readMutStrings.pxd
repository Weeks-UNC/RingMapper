

cdef struct READ:
    int start
    int stop
    char* read
    char* muts
    char* qual


cdef READ parseLine(char*, int)


cdef void fillReadMut(int[:], int[:], READ, int, int)

