###  TEMPLATE.txt.tpl; coding: utf-8 ---

#  Author(s): Christophe Prud'homme <christophe.prudhomme@feelpp.org>
#       Date: 2020-06-03
#
#  Copyright (C) 2020 Cemosis
#
# Distributed under the GPL(GNU Public License):
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#
feelpp_add_application(wave SRCS wave.cpp TESTS INSTALL )
feelpp_add_application(elasticity SRCS elasticity.cpp TESTS INSTALL )
# feelpp_add_application(elasticity_cantilever SRCS elasticity_cantilever.cpp TESTS INSTALL )
feelpp_add_testcase(cases PREFIX fp CATEGORY feelpp)

# if(FEELPP_TOOLBOXES_FOUND)
#     feelpp_add_application(toolbox SRCS toolbox.cpp LINK_LIBRARIES Feelpp::feelpp_toolbox_electric_lib TESTS)
#     feelpp_add_testcase(cases)
# endif()

if ( SKBUILD_PROJECT_NAME )
    python_add_library(_wave MODULE _wave.cpp WITH_SOABI)
    target_link_libraries(_wave PRIVATE pybind11::headers Feelpp::feelpp)
    target_compile_definitions(_wave PRIVATE VERSION_INFO=${PROJECT_VERSION})

    install(TARGETS _wave DESTINATION 2023_m2_feelpp_sonorhc)
else()
    feelpp_add_pymodule(wave SRCS _wave.cpp DESTINATION ${PROJECT_NAME})
    file(GLOB PYFILES "${CMAKE_CURRENT_SOURCE_DIR}/2023_m2_feelpp_sonorhc/*.py")
    install(FILES ${PYFILES} DESTINATION ${FEELPP_PYTHON_MODULE_PATH}/${PROJECT_NAME})
endif()

