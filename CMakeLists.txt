

cmake_minimum_required(VERSION 3.27)
project(CC1PiProject CXX)

# Establecer la versión mínima de macOS

# Establecer la versión mínima de macOS
set(CMAKE_OSX_DEPLOYMENT_TARGET "13.5" CACHE STRING "Minimum macOS version")

# Configurar las opciones de compilación y enlace
set(CMAKE_CXX_STANDARD 11)  # Establecer el estándar de C++
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=13.5")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mmacosx-version-min=13.5")


#find_package(ROOT REQUIRED COMPONENTS Hist)
find_package(ROOT REQUIRED COMPONENTS Core RIO Tree Hist MathCore Physics)

add_executable(CC1PiProject
        main.cpp
        TreeManagementCode/read_tree.C
        TreeManagementCode/make_class_tree.C
        TrueStudies/Double_diff_xsec.C
        TrueStudies/Resonance_plot.C
        RecoStudies/TrackPurCom/TL_impact_on_number_events.C
        TrueStudies/Interaction_final_state_plot.C
        RecoStudies/MissRec/MissReconstructionClassifier_Reorganized.C
        RecoStudies/MissRec/Energy_comp.C
        Classes/SpacePoint.cpp
        Classes/SpacePoint.h
        Classes/Hit.cpp
        Classes/Hit.h
        Classes/Slice.cpp
        Classes/Slice.h
        Classes/CaloPoint.cpp
        Classes/CaloPoint.h
        RecoStudies/TrackPurCom/TL_impact_on_energy_spectrum.C
        RecoStudies/MissRec/vtx_refinement/display_refinement.C
        Classes/RecoAlgorithms/TPCPlot.cpp
        Classes/RecoAlgorithms/TPCPlot.h
        Classes/RecoAlgorithms/RefinedSlice.cpp
        Classes/RecoAlgorithms/RefinedSlice.h
        RecoStudies/MissRec/2DReco/pre_processing.cpp
        Classes/RecoAlgorithms/Cluster.cpp
        Classes/RecoAlgorithms/Cluster.h
        Classes/RecoAlgorithms/HoughAlgorithm.cpp
        Classes/RecoAlgorithms/HoughAlgorithm.h
        Classes/RecoAlgorithms/SPoint.cpp
        Classes/RecoAlgorithms/SPoint.h
        Classes/RecoAlgorithms/LineEquation.cpp
        Classes/RecoAlgorithms/LineEquation.h
        Classes/RecoAlgorithms/HoughLine.cpp
        Classes/RecoAlgorithms/HoughLine.h
        Classes/RecoAlgorithms/DBSCAN.cpp
        Classes/RecoAlgorithms/DBSCAN.h
        RecoStudies/FinalBackground/remaining_background_p_energy.C
        RecoStudies/event_display.C
        RecoStudies/Chi2_stuff/chi2_to_tl.C
        RecoStudies/Chi2_stuff/chi2_to_tl.C
        toy.cpp
)

target_link_libraries(CC1PiProject PUBLIC  ROOT::Core ROOT::RIO ROOT::Tree ROOT::Hist ROOT::MathCore ROOT::Physics)