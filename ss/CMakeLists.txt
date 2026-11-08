if(NOT TARGET integration)
    ariadne_add_submodule(integration)
endif()
if(NOT TARGET equation_solving)
    ariadne_add_submodule(equation-solving)
endif()
if(NOT TARGET system)
    ariadne_add_submodule(system)
endif()

# The testing framework is a dependency of test executables only.
if(NOT TARGET ariadne-testing)
    ariadne_add_submodule(testing)
endif()

# Explicit development-only dependency for benchmark CLI composition.
if(NOT TARGET threading)
    ariadne_add_submodule(threading)
endif()
