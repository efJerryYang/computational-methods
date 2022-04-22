#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "GR::GKS" for configuration "Release"
set_property(TARGET GR::GKS APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GR::GKS PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libGKS.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libGKS.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS GR::GKS )
list(APPEND _IMPORT_CHECK_FILES_FOR_GR::GKS "${_IMPORT_PREFIX}/lib/libGKS.dll.a" "${_IMPORT_PREFIX}/bin/libGKS.dll" )

# Import target "GR::GR" for configuration "Release"
set_property(TARGET GR::GR APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GR::GR PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libGR.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libGR.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS GR::GR )
list(APPEND _IMPORT_CHECK_FILES_FOR_GR::GR "${_IMPORT_PREFIX}/lib/libGR.dll.a" "${_IMPORT_PREFIX}/bin/libGR.dll" )

# Import target "GR::GR3" for configuration "Release"
set_property(TARGET GR::GR3 APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GR::GR3 PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libGR3.dll.a"
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "GR::GR"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libGR3.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS GR::GR3 )
list(APPEND _IMPORT_CHECK_FILES_FOR_GR::GR3 "${_IMPORT_PREFIX}/lib/libGR3.dll.a" "${_IMPORT_PREFIX}/bin/libGR3.dll" )

# Import target "GR::GRM" for configuration "Release"
set_property(TARGET GR::GRM APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GR::GRM PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libGRM.dll.a"
  IMPORTED_LINK_DEPENDENT_LIBRARIES_RELEASE "GR::GR;GR::GR3"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libGRM.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS GR::GRM )
list(APPEND _IMPORT_CHECK_FILES_FOR_GR::GRM "${_IMPORT_PREFIX}/lib/libGRM.dll.a" "${_IMPORT_PREFIX}/bin/libGRM.dll" )

# Import target "GR::qt5gr" for configuration "Release"
set_property(TARGET GR::qt5gr APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(GR::qt5gr PROPERTIES
  IMPORTED_IMPLIB_RELEASE "${_IMPORT_PREFIX}/lib/libqt5gr.dll.a"
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/bin/libqt5gr.dll"
  )

list(APPEND _IMPORT_CHECK_TARGETS GR::qt5gr )
list(APPEND _IMPORT_CHECK_FILES_FOR_GR::qt5gr "${_IMPORT_PREFIX}/lib/libqt5gr.dll.a" "${_IMPORT_PREFIX}/bin/libqt5gr.dll" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
