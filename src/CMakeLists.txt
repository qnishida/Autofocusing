# 実行ファイルを作成
add_executable(${TARGET_NAME1} cal_ccf.cpp station_info.cpp util.cpp calTT.cpp)
add_executable(${TARGET_NAME2} cal_ccf_eq.cpp station_info.cpp util.cpp calTT.cpp)

foreach(target ${TARGET_NAME1} ${TARGET_NAME2})
  target_include_directories(${target} PRIVATE ${DELTAP_INCLUDE_DIRS})
  target_compile_options(${target} PRIVATE ${HDF5_DEFINITIONS} ${DELTAP_COMPILE_OPTIONS})
  target_link_libraries(${target} PRIVATE ${DELTAP_LIBS})
endforeach()
