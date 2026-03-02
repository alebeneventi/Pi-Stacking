# =========================
# Quantum Molecular Docking
# =========================

# Configuration
ENV_NAME := qmol_env
PYTHON_VERSION := 3.13.2
SRC_DIR := src
REQUIREMENTS := requirements.txt

# Default target
.DEFAULT_GOAL := help

# Environment Setup
# -----------------

create: ## Create conda environment with Python
	@echo "Creating conda environment: $(ENV_NAME)"
	conda create -y -n $(ENV_NAME) python=$(PYTHON_VERSION)
	@echo "Environment '$(ENV_NAME)' created successfully"

install: ## Install Python dependencies from requirements.txt
	@echo "Installing dependencies in $(ENV_NAME)"
	conda run -n $(ENV_NAME) pip install -r $(REQUIREMENTS)
	@echo "Dependencies installed successfully"

setup: create install ## Complete environment setup (create + install)

# Data Processing
# ---------------

preprocess: ## Run preprocessing on a single protein-ligand pair
	@echo "Running single-pair preprocessing"
	conda run -n $(ENV_NAME) python3 $(SRC_DIR)/preprocess.py

# Molecular Docking - Simulated Annealing
# ----------------------------------------

run-sa: ## Run simulated annealing on single protein-ligand pair
	@echo "Running single-pair simulated annealing"
	conda run -n $(ENV_NAME) python3 $(SRC_DIR)/dock_sa.py


# Environment Management
# ----------------------

env-info: ## Display environment information
	@echo "Environment: $(ENV_NAME)"
	@echo "Python Version: $(PYTHON_VERSION)"
	@conda info --envs | grep $(ENV_NAME) || echo "Environment not found"

env-export: ## Export environment to environment.yml
	conda env export -n $(ENV_NAME) > environment.yml
	@echo "Environment exported to environment.yml"

# Cleanup
# -------

clean: ## Remove Python cache files and temporary data
	@echo "Cleaning Python cache files"
	find . -type d -name "__pycache__" -delete
	@echo "Cache files cleaned"

clean-env: ## Remove conda environment
	@echo "Removing conda environment: $(ENV_NAME)"
	conda remove --name $(ENV_NAME) --all
	@echo "Environment removed"

clean-all: clean clean-env ## Full cleanup (cache + environment)

# Utilities
# ---------

help: ## Display this help message
	@echo "Quantum Molecular Docking - Available Commands:"
	@echo "================================================================================================="
	@awk 'BEGIN {FS = ":.*##"} /^[a-zA-Z_-]+:.*##/ { printf "  \033[36m%-26s\033[0m %s\n", $$1, $$2 }' $(MAKEFILE_LIST)
	@echo ""
	@echo "Quick Start:"
	@echo "  make setup                  # Create environment and install dependencies"
	@echo "  make preprocess             # Prepare data"
	@echo "  make run-sa                 # Run simulated annealing"
	@echo "  make clean-all              # Complete cleanup"
	@echo "================================================================================================="

# Phony targets
.PHONY: create install setup preprocess run-sa env-info env-export clean clean-env clean-all help