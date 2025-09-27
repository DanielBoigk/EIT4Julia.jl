module HilbertSpaces

export AbstractHilbertSpace,AbstractGalerkinSolver, AbstractBoundaryPair


abstract type AbstractHilbertSpace end # This is supposed to hold all information about the space:
abstract type AbstractGalerkinSolver end # This is supposed to hold all the information that the solver needs
abstract type AbstractBoundaryPair end # This is supposed to halod all the information for a voltage-current boundary pair.



end # Module
    
