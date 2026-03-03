# Check if Anaconda is installed
if ! command -v conda &> /dev/null; then
    echo "Anaconda is not installed. Please install Anaconda and try again."
    exit 1
fi

# Check if the environment is already installed
if [ ! -d $(conda info --base)/envs/raisd-ai ]; then
        conda env create -f raisd-ai.yml
fi

# download and compile RAiSD-AI and add it to the PATH
if [ ! -d ~/software/RAiSD-AI ]; then
    eval "$(conda shell.bash hook)"
    conda activate raisd-ai
    mydir=$(pwd)
    cd ~
    mkdir -p software
    cd software
    mkdir RAiSD-AI
    cd RAiSD-AI
    wget https://github.com/alachins/RAiSD-AI/archive/refs/heads/master.zip
    unzip master.zip && rm master.zip
    mv RAiSD-AI-master/* . && rm -r RAiSD-AI-master

    source ./compile-RAiSD-AI.sh
    if ! grep -q 'export PATH=.*~/software/RAiSD-AI' ~/.bashrc; then
        echo 'export PATH=$PATH:~/software/RAiSD-AI' >> ~/.bashrc
    fi
    if ! grep -q 'export PATH=.*~/software/RAiSD-AI' ~/.zshrc; then
        echo 'export PATH=$PATH:~/software/RAiSD-AI' >> ~/.zshrc
    fi

    source ~/.bashrc
    source ~/.zshrc
    cd $mydir
fi
