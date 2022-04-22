EXTRA_INCLUDEPATH += /workspace/srcdir/qtbase-everywhere-src-5.15.2/include/QtANGLE
host_build {
    QT_CPU_FEATURES.x86_64 = mmx sse sse2
} else {
    QT_CPU_FEATURES.x86_64 = mmx sse sse2
}
QT.global_private.enabled_features = sse2 alloca_malloc_h alloca dbus gui network relocatable sql system-zlib testlib widgets xml
QT.global_private.disabled_features = alloca_h android-style-assets avx2 private_tests dbus-linked dlopen gc_binaries intelcet libudev posix_fallocate reduce_exports reduce_relocations release_tools stack-protector-strong zstd
QT_COORD_TYPE = double
QMAKE_LIBS_ZLIB = -lz
CONFIG += cross_compile sse2 aesni sse3 ssse3 sse4_1 sse4_2 compile_examples largefile precompile_header rdrnd shani x86SimdAlways
QT_BUILD_PARTS += libs
QT_HOST_CFLAGS_DBUS += 
