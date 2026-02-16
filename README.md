Need to pull and compile n2dataread. This creates a static & dynamic library that this software needs to link with.
    I had to make some changes for mac compatability for n2dataread to compile
    
        CMakeLists.txt
            +set(CMAKE_C_FLAGS "-g -O0 -Wall -Wextra -pedantic -std=gnu1x -fdiagnostics-color=always -Wno-unused-function -I/usr/local/include -L/usr/local/lib -I/opt/homebrew/include -L/opt/homebrew/lib")
            +target_include_directories(N2readData_static PUBLIC ${CMAKE_SOURCE_DIR})
            +target_include_directories(N2readData_dyn    PUBLIC ${CMAKE_SOURCE_DIR})

            COMMENTED OUT ALL THE STUFF I DON'T NEED (AMI,EXECUTABLES)

            +install(TARGETS N2readData_dyn N2readData_static
                   LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/../install
                   ARCHIVE DESTINATION ${CMAKE_SOURCE_DIR}/../install
                   )

        src/ListHash/ListHash.h:
            +#ifdef __APPLE__
            +typedef int (*__compar_fn_t)(const void*, const void*);
            +#endif

        src/SimpleLog/SimpleLog.h
            ...
             #define SL_MSG_MAX 1024        // Max size of messages written to SimpleLog_Write
            +#ifdef __cplusplus^M
            +extern "C" {^M
            +#endif^M
             extern void SimpleLog_Setup(const char *PathName,
                                                                    const char *TimeFormat,
                                                                    const int NoRepeatLastN,
             ...
             extern void SimpleLog_Write(const int Level,
             extern void SimpleLog_Flush(void);
             extern void SimpleLog_Free(void);

            +#ifdef __cplusplus^M
            +}^M
            +#endif^M
            +^M
             #if __STDC_VERSION__ >= 199901L
                    // Even simpler: This is a possible shortcut in C99 only
                    #define SLOG(Level, ...) SimpleLog_Write(Level, SL_ORIGIN, __VA_ARGS__)
             ...


Also requires ROOT

I assume it is up one directory for my CMake configuration:

some_directory/
    n2dataread/
    data_reduction/CMakeLists.txt

