macro(set_abs_paths DEST)
	set(FILES "${ARGV}")
	list(REMOVE_AT FILES 0) # get rid of dest value
	foreach(_FILE ${FILES})
		get_filename_component(_TMP ${_FILE} ABSOLUTE)
		list(APPEND ${DEST} ${_TMP})
	endforeach(_FILE ${FILES})
endmacro(set_abs_paths DEST)
