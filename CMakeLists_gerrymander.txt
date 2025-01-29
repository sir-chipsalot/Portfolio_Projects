cmake_minimum_required( VERSION 3.20 )
project( gerrymander VERSION 1.0 )

message( "Using sources: gerrymander.cpp" )
add_executable( gerrymander )
target_sources( gerrymander PRIVATE gerrymander.cpp )

target_compile_features( gerrymander PRIVATE cxx_std_17 )
install( TARGETS gerrymander DESTINATION . )
