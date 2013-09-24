# Build specifications for the contents of the hardware package
#
# (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>

#
# Collect all code in this part of the source tree
#

add_library(hardware
	system.cpp
	device.cpp
	device_info.cpp
	profiling_data.cpp
	synchronization_event.cpp
)
target_link_libraries(hardware optimal)
add_modules(hardware
	buffers
	code
	transfer
)

#
# Definition of tets
#

add_unit_test(hardware::System system_test)
add_unit_test(hardware::Device device_test)
add_unit_test(hardware::ProfilingData profiling_data_test)
add_unit_test(hardware::SynchronizationEvent synchronization_event_test)
