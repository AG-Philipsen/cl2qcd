find_package( Nettle REQUIRED )
include_directories(SYSTEM ${Nettle_INCLUDE_DIRS})

# Build one crypto lib
add_library(crypto
	md5.cpp
)
target_link_libraries(crypto
	${Nettle_LIBRARIES}
)

add_executable(cryptoTest test.cpp)
target_link_libraries(cryptoTest ${Boost_LIBRARIES} crypto)
add_test(Crypto cryptoTest)
