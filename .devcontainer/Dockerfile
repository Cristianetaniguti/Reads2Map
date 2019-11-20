FROM python:3.7.3

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Or your actual UID, GID on Linux if not the default 1000
ARG USERNAME=vscode
ARG USER_UID=1000
ARG USER_GID=$USER_UID

# Configure apt and install packages
RUN apt-get update \
    && apt-get -y install --no-install-recommends apt-utils dialog 2>&1 \
    && apt-get -y install git procps lsb-release bash-completion shellcheck \
    && pip --disable-pip-version-check --no-cache-dir install pylint flake8 snoop \
    && groupadd --gid $USER_GID $USERNAME \
    && useradd -s /bin/bash --uid $USER_UID --gid $USER_GID -m $USERNAME \
    && apt-get install -y sudo \
    && echo $USERNAME ALL=\(root\) NOPASSWD:ALL > /etc/sudoers.d/$USERNAME \
    && chmod 0440 /etc/sudoers.d/$USERNAME \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*


ENV DOCKERIZE_VERSION v0.6.1
RUN wget https://github.com/jwilder/dockerize/releases/download/$DOCKERIZE_VERSION/dockerize-linux-amd64-$DOCKERIZE_VERSION.tar.gz \
    && tar -C /usr/local/bin -xzvf dockerize-linux-amd64-$DOCKERIZE_VERSION.tar.gz \
    && rm dockerize-linux-amd64-$DOCKERIZE_VERSION.tar.gz

# RUN apk add --no-cache tzdata
ENV TZ America/Sao_Paulo

# git prompt
RUN git clone https://github.com/magicmonty/bash-git-prompt.git /home/vscode/.bash-git-prompt --depth=1
RUN echo "\nif [ -f '/home/vscode/.bash-git-prompt/gitprompt.sh' ]; then\nGIT_PROMPT_ONLY_IN_REPO=1\nsource /home/vscode/.bash-git-prompt/gitprompt.sh\nfi" >> /home/vscode/.bashrc


# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=
