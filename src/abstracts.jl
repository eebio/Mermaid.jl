using DocStringExtensions

"""
$(TYPEDEF)
"""
abstract type AbstractComponent end

"""
$(TYPEDEF)
"""
abstract type AbstractTimeIndependentComponent <: AbstractComponent end

"""
$(TYPEDEF)
"""
abstract type AbstractTimeDependentComponent <: AbstractComponent end

"""
$(TYPEDEF)
"""
abstract type AbstractComponentIntegrator end

"""
$(TYPEDEF)
"""
abstract type AbstractMermaidSolver end

"""
$(TYPEDEF)
"""
abstract type AbstractMermaidIntegrator end

"""
$(TYPEDEF)
"""
abstract type AbstractMermaidProblem end

"""
$(TYPEDEF)
"""
abstract type AbstractMermaidSolution end

"""
$(TYPEDEF)
"""
abstract type AbstractConnectedVariable end

"""
$(TYPEDEF)
"""
abstract type AbstractConnector end
