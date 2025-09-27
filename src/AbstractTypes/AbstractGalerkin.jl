module HilbertSpaces

export AbstractHilbertSpace,AbstractGalerkinSolver, AbstractBoundaryPair


abstract type AbstractHilbertSpace end # This is supposed to hold all information about the space
abstract type AbstractGalerkinSolver end # This is supposed to hold all the information that the solver needs
abstract type AbstractAdjointSolver end # This is supposed to be all needed to get a gradient δσ for the update of σ 
abstract type AbstractBoundaryPair end # This is supposed to halod all the information for a voltage-current boundary pair.



end # Module
    
