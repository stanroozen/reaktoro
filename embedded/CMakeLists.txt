# Recursively collect runtime resource files only.
file(GLOB_RECURSE FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    databases/*.dat
    databases/*.json
    databases/*.txt
    databases/*.xlsm
    databases/*.xml
    databases/*.yaml
    interpolation/*.txt
    params/*.json
    params/*.yaml)

# Create a resource library containing embedded document files
cmrc_add_resource_library(ReaktoroEmbedded
    ALIAS Reaktoro::Embedded
    PREFIX embedded
    ${FILES}
)

# Set some target properties
set_target_properties(ReaktoroEmbedded PROPERTIES
    POSITION_INDEPENDENT_CODE ON)
