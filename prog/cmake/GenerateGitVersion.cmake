find_package(Git)

if(GIT_FOUND)
	execute_process( COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:%H
	                 WORKING_DIRECTORY ${SOURCE_DIR}
	                 OUTPUT_VARIABLE GIT_COMMIT_ID
	                 OUTPUT_STRIP_TRAILING_WHITESPACE )
	execute_process( COMMAND ${GIT_EXECUTABLE} diff --shortstat
	                 WORKING_DIRECTORY ${SOURCE_DIR}
	                 OUTPUT_VARIABLE GIT_DIRTY
	                 OUTPUT_STRIP_TRAILING_WHITESPACE )
	if( GIT_DIRTY )
		set( GIT_COMMIT_ID "${GIT_COMMIT_ID}*" )
	endif( GIT_DIRTY )
	message("Git commit id: ${GIT_COMMIT_ID}")
else(GIT_FOUND)
	message(WARNING "Failed to find git, can't write git commit id to executables")
	set(GIT_COMMIT_ID "N/A")
endif(GIT_FOUND)
file(WRITE gitcommitid.h.new "#define GIT_COMMIT_ID \"${GIT_COMMIT_ID}\"\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
		gitcommitid.h.new gitcommitid.h)
