#!/bin/bash
set -e

echo "Starting installation for SPORTA..."

# Function to check if a command exists
command_exists () {
    type "$1" &> /dev/null ;
}

# 1. System Dependency Check (Informational)
echo "Checking system dependencies..."
if ! command_exists pkg-config; then
    echo "Warning: pkg-config not found. Build might fail if libraries aren't in standard paths."
fi

# 2. Setup Python Environment
if [ -d "venv" ]; then
    echo "Virtual environment 'venv' already exists. Using it."
else
    echo "Creating virtual environment 'venv'..."
    python3 -m venv venv
fi

source venv/bin/activate

# 3. Upgrade pip and build tools
echo "Upgrading pip and setuptools..."
pip install --upgrade pip setuptools wheel

# 4. Install the package
echo "Installing SPORTA and dependencies..."

# Install the package (builds the C extension)
pip install .

echo "--------------------------------------------------------"
echo "Installation complete!"
echo "To run the application, simply define the environment and run 'sporta':"
echo "  source venv/bin/activate"
echo "  sporta"
echo ""
echo "Or directly:"
echo "  ./venv/bin/sporta"
echo "--------------------------------------------------------"
