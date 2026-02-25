#!/bin/bash

# FastJet Universal Installer Script
# Installs either 32-bit or 64-bit version of FastJet + FastJet Contrib
# Usage: ./fastjet_install.sh [x32|x64]
# If no argument provided, script will prompt for architecture choice

FASTJET_VERSION=3.4.2
FASTJET_CONTRIB_VERSION=1.046

# Function to display usage information
usage() {
    echo "Usage: $0 [x32|x64]"
    echo "  x32: Install 32-bit version (required for RCAS & PDSF)"
    echo "  x64: Install 64-bit version"
    echo "  If no argument provided, script will prompt for choice"
    exit 1
}

# Function to prompt user for architecture choice
prompt_architecture() {
    echo "Select architecture to install:"
    echo "1) x32 (32-bit) - Required for RCAS & SDCC"
    echo "2) x64 (64-bit)"
    read -p "Enter choice (1 or 2): " choice

    case $choice in
    1)
        ARCH="x32"
        ;;
    2)
        ARCH="x64"
        ;;
    *)
        echo "Invalid choice. Exiting."
        exit 1
        ;;
    esac
}

# Function to install FastJet with 32-bit configuration
install_fastjet_x32() {
    echo "Installing FastJet $FASTJET_VERSION (32-bit)..."

    # Configure with 32-bit flags
    ./configure --prefix=$PWD/../fastjet-install \
        CXXFLAGS="-m32 -fPIC -fno-inline" \
        CFLAGS="-m32 -fPIC -fno-inline" \
        LDFLAGS="-m32"
}

# Function to install FastJet with 64-bit configuration
install_fastjet_x64() {
    echo "Installing FastJet $FASTJET_VERSION (64-bit)..."

    # Setup 64-bit environment if needed
    if command -v "setup" &>/dev/null; then
        setup 64b
    fi

    # Configure without special flags
    ./configure --prefix=$PWD/../fastjet-install
}

# Function to install FastJet Contrib with 32-bit configuration
install_contrib_x32() {
    echo "Installing FastJet Contrib $FASTJET_CONTRIB_VERSION (32-bit)..."

    ./configure --fastjet-config=$PWD/../fastjet-$FASTJET_VERSION/fastjet-config \
        --prefix=$PWD/../fastjet-install \
        CXXFLAGS="-m32 -fPIC -fno-inline" \
        CFLAGS="-m32 -fPIC -fno-inline" \
        LDFLAGS="-m32"
}

# Function to install FastJet Contrib with 64-bit configuration
install_contrib_x64() {
    echo "Installing FastJet Contrib $FASTJET_CONTRIB_VERSION (64-bit)..."

    ./configure --fastjet-config=$PWD/../fastjet-$FASTJET_VERSION/fastjet-config \
        --prefix=$PWD/../fastjet-install
}

# Function to setup environment variables
setup_environment() {
    echo "Setting up environment variables..."

    # Export variables for current session
    export FASTJET=$install_dir/fastjet-install
    export FASTJET_CONTRIB=$install_dir/fjcontrib-$FASTJET_CONTRIB_VERSION

    # Add to .bashrc for persistent setup
    echo "export FASTJET='$install_dir/fastjet-install'" >>~/.bashrc
    echo "export FASTJET_CONTRIB='$install_dir/fjcontrib-$FASTJET_CONTRIB_VERSION'" >>~/.bashrc

    echo "Environment variables added to ~/.bashrc"
}

# Function to create symbolic links
create_symlinks() {
    echo "Creating symbolic links for include paths..."
    cd $current_dir

    if [ -d "$FASTJET/include/siscone" ]; then
        ln -sf $FASTJET/include/siscone siscone
        echo "Created symlink: siscone -> $FASTJET/include/siscone"
    fi

    if [ -d "$FASTJET/include/fastjet" ]; then
        ln -sf $FASTJET/include/fastjet fastjet
        echo "Created symlink: fastjet -> $FASTJET/include/fastjet"
    fi
}

# Main installation function
main() {
    # Parse command line arguments
    if [ $# -eq 0 ]; then
        prompt_architecture
    elif [ $# -eq 1 ]; then
        case $1 in
        x32 | 32)
            ARCH="x32"
            ;;
        x64 | 64)
            ARCH="x64"
            ;;
        -h | --help)
            usage
            ;;
        *)
            echo "Error: Invalid argument '$1'"
            usage
            ;;
        esac
    else
        echo "Error: Too many arguments"
        usage
    fi

    # Set installation directory based on architecture
    install_dir="/gpfs01/star/pwg/$USER/install/$ARCH"
    current_dir=$(pwd)

    echo "Starting FastJet installation for $ARCH architecture..."
    echo "Installation directory: $install_dir"

    # Create installation directory
    mkdir -p $install_dir
    cd $install_dir

    # Download and extract FastJet
    echo "Downloading FastJet $FASTJET_VERSION..."
    curl -O https://fastjet.fr/repo/fastjet-$FASTJET_VERSION.tar.gz
    tar zxvf fastjet-$FASTJET_VERSION.tar.gz
    cd fastjet-$FASTJET_VERSION/

    # Install FastJet based on architecture
    if [ "$ARCH" = "x32" ]; then
        install_fastjet_x32
    else
        install_fastjet_x64
    fi

    # Build and install FastJet
    make
    make check
    make install

    # Download and extract FastJet Contrib
    cd ..
    echo "Downloading FastJet Contrib $FASTJET_CONTRIB_VERSION..."
    wget http://fastjet.hepforge.org/contrib/downloads/fjcontrib-$FASTJET_CONTRIB_VERSION.tar.gz
    tar zxvf fjcontrib-$FASTJET_CONTRIB_VERSION.tar.gz
    cd fjcontrib-$FASTJET_CONTRIB_VERSION

    # Install FastJet Contrib based on architecture
    if [ "$ARCH" = "x32" ]; then
        install_contrib_x32
    else
        install_contrib_x64
    fi

    # Build and install FastJet Contrib
    make
    make install

    # Additional steps for 32-bit installation
    if [ "$ARCH" = "x32" ]; then
        echo "Building fragile shared libraries for 32-bit..."
        make fragile-shared
        make fragile-shared-install

        # Unset compiler flags
        unset CXXFLAGS CFLAGS LDFLAGS
        echo "Unset 32-bit compiler flags"
    fi

    # Setup environment variables
    setup_environment

    # Source .bashrc to update current session
    source ~/.bashrc

    # Create symbolic links
    create_symlinks

    echo ""
    echo "FastJet installation completed successfully!"
    echo "Architecture: $ARCH"
    echo "FastJet version: $FASTJET_VERSION"
    echo "FastJet Contrib version: $FASTJET_CONTRIB_VERSION"
    echo "Installation path: $install_dir/fastjet-install"
    echo ""
    echo "Environment variables have been added to ~/.bashrc"
    echo "Symbolic links created in: $current_dir"
    echo ""
    echo "You may need to restart your shell or run 'source ~/.bashrc' to use FastJet"
}

# Run main function with all arguments
main "$@"
