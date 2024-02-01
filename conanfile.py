from conans import ConanFile, CMake, tools


class SmsppConan(ConanFile):
    name = "milpsolver"
    version = "0.7.1"
    description = "A MILP Solver for SMS++"
    topics = ("conan", "smspp", "milpsolver")
    url = "https://gitlab.com/smspp/milpsolver"
    homepage = "https://gitlab.com/smspp/milpsolver"
    license = "GPL-3.0-only"
    generators = "cmake"

    settings = "os", "arch", "compiler", "build_type"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {"shared": False, "fPIC": True}

    requires = (
        "smspp/0.5.2@smspp/testing",
        "CPLEX/12.10"
    )

    exports_sources = [
        "CMakeLists.txt",
        "src/*",
        "include/*",
        "tools/*",
        "cmake/*"
    ]

    def source(self):
        tools.replace_in_file(
            "CMakeLists.txt",
            '''LANGUAGES C CXX)''',
            '''LANGUAGES C CXX)\n''' +
            '''include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)\n''' +
            '''conan_basic_setup()'''
        )

    def _configure_cmake(self):
        cmake = CMake(self)
        cmake.definitions["BUILD_TESTING"] = False
        # SCIP not supported in Conan, yet
        cmake.definitions["MILPSolver_USE_SCIP"] = False
        cmake.configure()
        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()

    def package(self):
        self.copy(pattern="LICENSE", dst="licenses")
        cmake = self._configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.includedirs = ["include", "include/SMS++"]
        self.cpp_info.libs = ["MILPSolver"]
