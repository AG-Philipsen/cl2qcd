include(FindGit REQUIRED)

execute_process( COMMAND ${GIT_EXECUTABLE} log -n 1 --pretty=format:%H
                 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                 OUTPUT_VARIABLE GIT_COMMIT_ID
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
execute_process( COMMAND ${GIT_EXECUTABLE} diff --shortstat
                 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
                 OUTPUT_VARIABLE GIT_DIRTY
                 OUTPUT_STRIP_TRAILING_WHITESPACE )
if( GIT_DIRTY )
	set( GIT_COMMIT_ID "${GIT_COMMIT_ID}*" )
endif( GIT_DIRTY )

file(WRITE gitcommitid.h.new "#define GIT_COMMIT_ID \"${GIT_COMMIT_ID}\"\n")
execute_process(COMMAND ${CMAKE_COMMAND} -E copy_if_different
		gitcommitid.h.new gitcommitid.h)
