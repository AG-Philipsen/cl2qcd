# MD5 is endian-dependent
include(TestBigEndian REQUIRED)
test_big_endian(BIG_ENDIAN)
if(BIG_ENDIAN)
	add_definitions( "-DMD5_BIG_ENDIAN=true" )
endif(BIG_ENDIAN)

# Define sources for MD5
set(MD5_SRCS
	md5.c
)

# Build one crypto lib
add_library(crypto ${MD5_SRCS})
