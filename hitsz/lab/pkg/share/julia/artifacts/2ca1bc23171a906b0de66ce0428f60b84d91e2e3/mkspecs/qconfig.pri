host_build {
    QT_ARCH = x86_64
    QT_BUILDABI = x86_64-little_endian-lp64
    QT_TARGET_ARCH = x86_64
    QT_TARGET_BUILDABI = x86_64-little_endian-lp64
} else {
    QT_ARCH = x86_64
    QT_BUILDABI = x86_64-little_endian-lp64
}
QT.global.enabled_features = shared cross_compile shared c++11 c++14 c++17 c++1z c99 c11 thread future concurrent signaling_nan
QT.global.disabled_features = framework rpath appstore-compliant debug_and_release simulator_and_device build_all c++2a pkg-config force_asserts separate_debug_info static
QT_CONFIG += shared shared release c++11 c++14 c++17 c++1z concurrent dbus no-pkg-config stl
CONFIG += shared cross_compile shared release no_plugin_manifest
QT_VERSION = 5.15.2
QT_MAJOR_VERSION = 5
QT_MINOR_VERSION = 15
QT_PATCH_VERSION = 2
QT_GCC_MAJOR_VERSION = 8
QT_GCC_MINOR_VERSION = 1
QT_GCC_PATCH_VERSION = 0
QT_EDITION = OpenSource
