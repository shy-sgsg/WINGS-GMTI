cmake_minimum_required(VERSION 3.18)

# ★ 关键修改 1: 必须显式启用 CUDA 语言支持，否则 nvcc 不会工作
project(GMTI_core LANGUAGES CXX CUDA)

# ===== 编译选项 =====
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ★ 关键修改 2: 开启 CUDA 分离编译
# 这样 GMTIProcessor 类的成员函数才能跨 .cpp 和 .cu 文件链接
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)

# 设置 CUDA 架构（根据你的显卡调整，75对应RTX20系列, 80对应A100, 86对应RTX30系列）
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 75 CACHE STRING "CUDA architectures")
endif()

set(LIBOMP_ROOT "/opt/homebrew/opt/libomp")  # OpenMP for Apple M1/M2

if(APPLE)
  set(CMAKE_OSX_ARCHITECTURES "arm64" CACHE STRING "" FORCE)
endif()

if(MSVC)
  add_compile_options(/W4 /permissive-)
else()
  add_compile_options(
    $<$<COMPILE_LANGUAGE:CXX>:-Wall>
    $<$<COMPILE_LANGUAGE:CXX>:-Wextra>
    $<$<COMPILE_LANGUAGE:CXX>:-Wpedantic>
  )
endif()

# ★ 开关选项
option(ENABLE_OPENMP "Enable OpenMP parallelization" OFF)
option(ENABLE_DEBUG_LOG "Force define DEBUG macro for verbose logging" OFF)

# ===== 依赖库查找 =====

# 1) CUDA Toolkit (现代方式查找)
find_package(CUDAToolkit REQUIRED)

# 3) OpenMP
if(ENABLE_OPENMP)
  find_package(OpenMP)
endif()

# 4) FFTW3 (保留你原来的复杂查找逻辑)
set(FFTW3_TARGET "")
set(FFTW3_THREADS_TARGET "")

find_package(FFTW3 QUIET COMPONENTS FFTW3 FFTW3_THREADS)
if(FFTW3_FOUND)
  if(TARGET FFTW3::FFTW3)
    add_library(FFTW3::fftw3 ALIAS FFTW3::FFTW3)
  endif()
  if(TARGET FFTW3::FFTW3_THREADS)
    add_library(FFTW3::fftw3_threads ALIAS FFTW3::FFTW3_THREADS)
  endif()
endif()

if(TARGET FFTW3::fftw3)
  set(FFTW3_TARGET FFTW3::fftw3)
elseif(TARGET fftw3::fftw3)
  add_library(FFTW3::fftw3 ALIAS fftw3::fftw3)
  set(FFTW3_TARGET FFTW3::fftw3)
endif()

if(TARGET FFTW3::fftw3_threads)
  set(FFTW3_THREADS_TARGET FFTW3::fftw3_threads)
elseif(TARGET fftw3::fftw3_threads)
  add_library(FFTW3::fftw3_threads ALIAS fftw3::fftw3_threads)
  set(FFTW3_THREADS_TARGET FFTW3::fftw3_threads)
endif()

# 兜底查找
if(NOT FFTW3_TARGET)
  find_library(FFTW3_LIB NAMES fftw3 PATHS /usr/lib64 /usr/lib /usr/local/lib /opt/homebrew/lib)
  if(FFTW3_LIB)
    add_library(FFTW3::fftw3 INTERFACE IMPORTED)
    target_link_libraries(FFTW3::fftw3 INTERFACE "${FFTW3_LIB}")
    set(FFTW3_TARGET FFTW3::fftw3)
  endif()
endif()

# 5) 线程库
find_package(Threads REQUIRED)

# TinyXML2 (needed by trackModule tool)

# ===== 源码配置 =====
set(SRC
  src/main.cpp
  src/alignFFTAndDBS.cpp
  src/clutter_CSI.cpp
  src/computeDynamicSupportDomain.cpp
  src/detectAndLocating.cpp
  src/loadXML.cpp
  src/processOnePeriod.cpp
  src/processPeriodsParallel.cpp
  src/pulseCompression.cpp
  src/rg_correct.cpp
  src/rangeCompress_cu.cu
  src/rangeCompress_device.cu
  src/writeresults.cpp
  src/kalman_track.cpp
  src/trackModule.cpp
  src/extractUTC.cpp
  src/makeLists.cpp
  src/pesudoTarget_Gen.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
  # ★ 确保这里的路径和你实际的文件名一致
  src/gpu_kernels.cu
)

# 生成可执行文件
add_executable(GMTI_core ${SRC})

add_executable(compare_range src/tools/compare_range.cpp)
target_include_directories(compare_range PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(compare_range PRIVATE
  ${FFTW3_TARGET}
  Threads::Threads
  CUDA::cufft
  CUDA::cudart
)
if(FFTW3_THREADS_TARGET)
  target_link_libraries(compare_range PRIVATE ${FFTW3_THREADS_TARGET})
endif()
target_sources(compare_range PRIVATE
  src/alignFFTAndDBS.cpp
  src/clutter_CSI.cpp
  src/computeDynamicSupportDomain.cpp
  src/detectAndLocating.cpp
  src/loadXML.cpp
  src/processOnePeriod.cpp
  src/processPeriodsParallel.cpp
  src/pulseCompression.cpp
  src/rg_correct.cpp
  src/rangeCompress_cu.cu
  src/rangeCompress_device.cu
  src/writeresults.cpp
  src/kalman_track.cpp
  src/extractUTC.cpp
  src/makeLists.cpp
  src/pesudoTarget_Gen.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
  src/gpu_kernels.cu
)

add_executable(convert_old_echo_to_new_format
  src/tools/convert_old_echo_to_new_format.cpp
  src/pulseCompression.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
)
target_include_directories(convert_old_echo_to_new_format PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(convert_old_echo_to_new_format PRIVATE
  ${FFTW3_TARGET}
  Threads::Threads
)
if(FFTW3_THREADS_TARGET)
  target_link_libraries(convert_old_echo_to_new_format PRIVATE ${FFTW3_THREADS_TARGET})
endif()

# ===== 包含目录与链接 =====

target_include_directories(GMTI_core PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CUDAToolkit_INCLUDE_DIRS}
)

# ★ 关键修改 3: 链接 CUDA 相关库
target_link_libraries(GMTI_core
  PRIVATE
    ${FFTW3_TARGET}
    Threads::Threads
    CUDA::cufft      # cuFFT 库
    CUDA::cudart     # CUDA Runtime 库
)

if(FFTW3_THREADS_TARGET)
  target_link_libraries(GMTI_core PRIVATE ${FFTW3_THREADS_TARGET})
endif()

# NOTE: GMTI_track target removed - track association now integrated into GMTI_core
# trackModule.cpp is now compiled as part of GMTI_core and called from main.cpp

# OpenMP 特殊处理 (Apple 兼容)
if(ENABLE_OPENMP)
  if(APPLE AND OpenMP_CXX_FOUND)
    target_compile_options(GMTI_core PRIVATE -Xpreprocessor -fopenmp)
    target_link_libraries(GMTI_core PRIVATE "${LIBOMP_ROOT}/lib/libomp.dylib")
    target_include_directories(GMTI_core PRIVATE "${LIBOMP_ROOT}/include")
  elseif(OpenMP_CXX_FOUND)
    target_link_libraries(GMTI_core PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

# ===== 宏定义与优化 =====
if(ENABLE_DEBUG_LOG)
  target_compile_definitions(GMTI_core PRIVATE DEBUG)
else()
  target_compile_definitions(GMTI_core PRIVATE
    $<$<CONFIG:Debug>:DEBUG>
    $<$<CONFIG:Release>:NDEBUG>
    $<$<CONFIG:RelWithDebInfo>:NDEBUG>
  )
endif()

# CUDA 编译优化选项
target_compile_options(GMTI_core PUBLIC
  $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math -Xcompiler -fPIC>
)
target_compile_options(GMTI_core PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>)

if(NOT MSVC)
  target_compile_options(GMTI_core PRIVATE
    $<$<CONFIG:Debug>:-O0 -g3 -fno-omit-frame-pointer>
    $<$<CONFIG:Release>:-O3>
  )
endif()

# ===== GMTI Pipe 可执行文件（Phase 1 框架）=====
# 与 GMTI_core 共享核心代码，但添加 pipe 通信层
set(GMTI_PIPE_SRC
  src/pipe/main.cpp
  src/pipe/MainCtrl.cpp
  src/pipe/PipeRW.cpp
  src/alignFFTAndDBS.cpp
  src/clutter_CSI.cpp
  src/computeDynamicSupportDomain.cpp
  src/detectAndLocating.cpp
  src/loadXML.cpp
  src/processOnePeriod.cpp
  src/processPeriodsParallel.cpp
  src/pulseCompression.cpp
  src/rg_correct.cpp
  src/rangeCompress_cu.cu
  src/rangeCompress_device.cu
  src/writeresults.cpp
  src/kalman_track.cpp
  src/extractUTC.cpp
  src/makeLists.cpp
  src/trackModule.cpp
  src/pesudoTarget_Gen.cpp
  src/tinystr.cpp
  src/tinyxml.cpp
  src/tinyxmlerror.cpp
  src/tinyxmlparser.cpp
  src/gpu_kernels.cu
)

add_executable(GMTI_pipe_core ${GMTI_PIPE_SRC})

target_include_directories(GMTI_pipe_core PRIVATE 
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/src/pipe
    ${CUDAToolkit_INCLUDE_DIRS}
)

target_link_libraries(GMTI_pipe_core
  PRIVATE
    ${FFTW3_TARGET}
    Threads::Threads
    CUDA::cufft
    CUDA::cudart
)

if(FFTW3_THREADS_TARGET)
  target_link_libraries(GMTI_pipe_core PRIVATE ${FFTW3_THREADS_TARGET})
endif()

# OpenMP 支持
if(ENABLE_OPENMP)
  if(APPLE AND OpenMP_CXX_FOUND)
    target_compile_options(GMTI_pipe_core PRIVATE -Xpreprocessor -fopenmp)
    target_link_libraries(GMTI_pipe_core PRIVATE "${LIBOMP_ROOT}/lib/libomp.dylib")
    target_include_directories(GMTI_pipe_core PRIVATE "${LIBOMP_ROOT}/include")
  elseif(OpenMP_CXX_FOUND)
    target_link_libraries(GMTI_pipe_core PRIVATE OpenMP::OpenMP_CXX)
  endif()
endif()

# 调试定义
if(ENABLE_DEBUG_LOG)
  target_compile_definitions(GMTI_pipe_core PRIVATE DEBUG)
else()
  target_compile_definitions(GMTI_pipe_core PRIVATE
    $<$<CONFIG:Debug>:DEBUG>
    $<$<CONFIG:Release>:NDEBUG>
    $<$<CONFIG:RelWithDebInfo>:NDEBUG>
  )
endif()

# CUDA 编译选项
target_compile_options(GMTI_pipe_core PUBLIC
  $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math -Xcompiler -fPIC>
)
target_compile_options(GMTI_pipe_core PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math>)

if(NOT MSVC)
  target_compile_options(GMTI_pipe_core PRIVATE
    $<$<CONFIG:Debug>:-O0 -g3 -fno-omit-frame-pointer>
    $<$<CONFIG:Release>:-O3>
  )
endif()