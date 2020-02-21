target_sources(mast
    PRIVATE
        ${CMAKE_CURRENT_LIST_DIR}/element_property_card_1D.cpp
        ${CMAKE_CURRENT_LIST_DIR}/element_property_card_1D.h
        ${CMAKE_CURRENT_LIST_DIR}/element_property_card_2D.cpp
        ${CMAKE_CURRENT_LIST_DIR}/element_property_card_2D.h
        ${CMAKE_CURRENT_LIST_DIR}/element_property_card_base.h
        ${CMAKE_CURRENT_LIST_DIR}/isotropic_element_property_card_3D.cpp
        ${CMAKE_CURRENT_LIST_DIR}/isotropic_element_property_card_3D.h
        ${CMAKE_CURRENT_LIST_DIR}/isotropic_material_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/isotropic_material_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/material_property_card_base.h
        ${CMAKE_CURRENT_LIST_DIR}/multilayer_1d_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/multilayer_1d_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/multilayer_2d_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/multilayer_2d_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/orthotropic_element_property_card_3D.cpp
        ${CMAKE_CURRENT_LIST_DIR}/orthotropic_element_property_card_3D.h
        ${CMAKE_CURRENT_LIST_DIR}/orthotropic_material_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/orthotropic_material_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_arbitrary_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_arbitrary_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_nparameter_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_nparameter_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_1parameter_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_1parameter_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_bar_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_bar_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_rod_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/solid_1d_rod_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/solid_2d_section_element_property_card.cpp
        ${CMAKE_CURRENT_LIST_DIR}/solid_2d_section_element_property_card.h
        ${CMAKE_CURRENT_LIST_DIR}/cross_section_property_pilkey.cpp
        ${CMAKE_CURRENT_LIST_DIR}/cross_section_property_pilkey.h)

# Install MAST headers for this directory.
install(DIRECTORY ./ DESTINATION include/property_cards
        FILES_MATCHING PATTERN "*.h")
