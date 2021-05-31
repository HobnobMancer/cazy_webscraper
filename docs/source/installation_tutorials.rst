================================================================
How to install ``cazy_webscraper``
================================================================

This page includes a step-by-step of how to check your computer setup meets the requirements to install ``cazy_webscraper``, and covers 3 different ways to install ``cazy_webscraper``. 

If you are experience using ``pip`` and ``bioconda`` to install Python packages, jump down to the relevant sections to get the command to install the packages.

If you are not experienced with installing programs via the command-line then read on!


Checking the requirements
****************************

To use ``cazy_webscraper`` you need to be using a POISx operating system (OS), Mac OS or a Linux emulator. This is because 
you need access to a Unix shell (a form of command-line) that writes in a language called Bash.

Linux
========
POISx includes operating systems such as Unix and it derivaties, Linux, and its derivatives such as Ubuntu. On Linux 
OS the default Unix shell is Bash.

Mac
=======
If you are working on a Mac then you already have a Unix shell program installed.

Mac OS Mojave or earlier releases
-------------------------------------

If you are running a Mac computer running Mac OS Mojave or earlier releases, the default Unix shell is using Bash. To access 
the Unix shell (or terminal) try one of the following:

* In Finder, select the Go menu, then select Utilities. Locate Terminal in the Utilities folder and open it
* Use the Mac ‘Spotlight’ computer search function. Search for: Terminal and press Return

To check if your machine is set up to use something other than Bash, type ``echo $SHELL`` in your terminal window.

Mac OS Catalina or later releases
-------------------------------------

For a Mac computer running macOS Catalina or later releases, or if you computer is set up to use something other 
than Bash, you can run Bash by opening a terminal, typing the command ``bash`` and hitting return.

Windows
===========

If you use Windows you can use a Linux emulator. Several emulators exist. Here we'll run through two different emulators you could use, 
and how to install them. These are not the only two and feel free to use which ever Linux emulator you wish!

Ubuntu
---------

Another Linux emulator that you can run on Windows is provided by Ubuntu. Ubuntu is a version of Linux with a graphical 
user interface (GUI), and they also provide a Windows emulator version that is fully supported, free and available via the `Microsoft app store <https://www.microsoft.com/en-gb/p/ubuntu-2004-lts/9n6svws3rx71#activetab=pivot:overviewtab>`_.

To install the Ubuntu emulator, navigate the `Microsoft store page <https://www.microsoft.com/en-gb/p/ubuntu-2004-lts/9n6svws3rx71#activetab=pivot:overviewtab>`_, and click install. 
This will handel all installation for you.

To start the terminal, open the Windows Start Menu and type 'Ubuntu', this will find the Ubuntu program. Click on Ubuntu and this will open a Ubuntu terminal. 

Ubuntu sets the Home as its own set of directories. To navigate to any hardrive in your system, open the Ubuntu terminal then type:  
``cd /mnt/<letter_of_harddrive (in lower case)>``, for example to access the C drive you would use ``cd /mnt/c``. The ``/mnt`` prefix 
accesses the Windows 'mounting' system, which is how it accesses harddrives within the computer.

Git for Windows
-----------------

One Linux emulator you can use if included in *Git for Windows*. The following installation instructions are adapted from 
Software Carpentry `lessons <https://carpentries.github.io/workshop-template/#shell>`_. You can find a video tutorial with a step-by-step guide `here <https://youtu.be/339AEqk9c-8>`_.

Installation instrucitions:
* Download the Git for Windows install from `here <https://gitforwindows.org/>`_.
* Run the installer
* Click the "Next" button four times (only two times if you already have Git installed). *You do not need to change anything the Information, location, components, and start menu screens.
* From the dropdown menu select "Use the Nano editor by default" (NOTE: you will need to scroll up to find it) and click "Next"
* On the page that says "Adjusting the name of the initial branch in new repositories", ensure that "Let Git decide" is selected.
* Ensure that "Git from the command line and also from 3rd-party software" is selected and click on "Next"
* Ensure that "Use the native Windows Secure Channel Library" is selected and click on "Next".
* Ensure that "Checkout Windows-style, commit Unix-style line endings" is selected and click on "Next".
* Ensure that "Use Windows' default console window" is selected and click on "Next".
* Ensure that "Default (fast-forward or merge) is selected and click "Next"
* Ensure that "Git Credential Manager Core" is selected and click on "Next".
* Ensure that "Enable file system caching" is selected and click on "Next".
* Click on "Install".
* Click on "Finish" or "Next".
* If your "HOME" environment variable is not set (or you don't know what this is):
* Open command prompt (Open the Windows Start Menu, then type 'cmd', and press Enter). Wait for the command promt window to open, and type the following line into the command prompt window exactly as shown:
``setx HOME "%USERPROFILE%"``. Press Enter, you should see ``SUCCESS: Specified value was saved``. Quit command prompt by typing `exit` then pressing Enter


The shell terminal
*********************
For more information on using the Unix shell checkout the content and lessons hosted at `SoftWare Carpentry <https://swcarpentry.github.io/shell-novice/01-intro/index.html>`_, which 
will walk through an introduction to the Unix shell and how the Unix shell can be used to optimise computational work.


Installing cazy_webscraper
******************************

There are three different ways to install ``cazy_webscraper``:  

* Via ``Bioconda``
* Via ``pip``
* From the source repository

Installing via ``Bioconda`` and ``pip` are by far the easiest methods of installation. Installing from source is not necessarily more difficult, but it takes a couple more steps.

Installing via ``Bioconda``
==============================

If you already have `conda` installed *and* the `bioconda` channel available, ``cazy_webscraper`` can be fully installed using the command:  

.. code-block:: bash

   conda install cazy_webscraper

If you do not have the `bioconda` channel available (you may discover this when trying to use the above command and the computer throws up the message ``PackagesNotFoundError: The following packages are not available from current channels``), you can install ``cazy_webscraper`` using the command:  

.. code-block:: bash

   conda install -c bioconda cazy_webscraper

If Conda is not installed, please see the Conda website for installation `instructions <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_.

.. warning::
   If you install ``cazy_webscraper`` via ``bioconda``, to invoke ``cazy_webscraper`` call it via the command-line using ``cazy_webscraper.py``.

Installing via ``pip``
==============================

``cazy_webscraper`` can also be installed via the `Python package index (Pypi)s <https://pypi.org/project/cazy-webscraper/>`_. To install ``cazy_webscraper`` use the command:  

.. code-block:: bash

   pip3 install cazy-webscraper

Now ``cazy_webscraper`` is fully installed and you can skip to the part of the tutorial that explains how to use it!

If ``pip`` is not installed, text will appear in the terminal telling you how to use `pip`. If `pip` is not install an error message will be displayed, stating that the computer could not find a command called `pip`. If this happens install pip using:  

.. code-block:: bash
   
   python get-pip.py

To update the version of `pip` install on your computer use:  


.. code-block:: bash
   
   python -m pip install --upgrade pip
 
 Then you should be able to install ``cazy_webscraper`` using the command provided above.
 
.. warning::
   If you install ``cazy_webscraper`` via ``pip``, to invoke ``cazy_webscraper`` call it via the command-line using ``cazy_webscraper.py``.

 
 
Installing from source: Installing using the files stored in the GitHub repository
====================================================================================


You can install ``cazy_webscraper`` directly from the GitHub repository. Open up your Bash terminal. Then we need to navigate to directory where you want to install ``cazy_webscraper``. To do this we will use the ``cd`` command.

Just like how the windows explorer points at a single directory at any time, and shows you the content of the directory, the terminal acts the same way. 
To check at what directory your terminal is pointed at (or looking at, type the command ``pwd`` and press enter. The terminal will then 
return the path of the directory at which it is currently looking at. For example, if the terminal is pointed at a directory called 'my_data' within another directory called 'Documents', located on the C drive, 
``pwd`` will return ``c/Documents/my_dir/``.

We can change directory using the 'change directory' command (``cd``). Continuing on from the example above, 
if we wanted to move from the directory 'my_dir' into the directory 'cazyme_research' located within it, and then into 
the directory 'cazy_dir' within that, we would use the command ``cd cazyme_research/cazy_dir``.

The change directory (``cd``) command is called then provided a path to the directory that we wish to 
have the terminal pointed at. The ``cd`` command starts at the directory the terminal is currently looking at, then 
follows the path we provide it. This is why to move from 'my_dir' > 'cazy_research' > 'cazy_dir' we can type 
``cd cazy_research/cazy_dir``, becuase the terminal will looking within the 'my_dir' directory for the 'cazy_research' directory.

Using the ``cd`` command navigate to the directory you wish to install ``cazy_webscraper``. 
**If the directory where you wish to install ``cazy_webscraper`` does not exist we can use the terminal to make it**. 
To use the terminal, first use the ``cd`` comamand to navigate to the parent directory of where you wish to house the 
directory that you will install ``cazy_webscraper``. Then call the 'make directory command' ``mkdir`` followed by the name 
you wish to give the directory. For example, once we have navigated to the 'cazy_dir', we can make the directory 
'cazyme_databases' by using ``mkdir cazyme_database``. We can then navigate into the 'cazyme_database' directory we have justed made 
by typing ``cd cazy_database`` into the terminal and hitting Return.

**Installing ``cazy_webscraper``**

First we clone the GitHub repository, by using the code:

.. code-block:: bash

   git clone https://github.com/HobnobMancer/cazy_webscraper 

This creates a new directory into the directory that the terminal is currently pointed at, called 
'cazy_webscraper'. The command also downloads all files in the GitHub repository, and writes them into 
the new 'cazy_webscraper' directory.

We then need to move into the 'cazy_webscraper' directory:

.. code-block:: bash

   cd cazy_webscraper

We then use the Python package manage ``pip`` to install ``cazy_webscraper``.

.. code-block:: bash

   pip3 install -e .

Do not forget the **-e** from this command, otherwise ``cazy_webscraper`` will not be installed correctly 
and you will run into constant issues when trying to use ``cazy_webscraper``.

**If you ever invoke ``cazy_webscraper`` and want to cancle the command, simple press the ``Ctrl`` and ``c`` keys at the same time.**

.. warning::
   If you install ``cazy_webscraper`` from source, to invoke ``cazy_webscraper`` you will need call Python3 followed by the path to the ``cazy_webscraper.py`` file. For example, if you are located in the root directory of the respository, you would use:  ``pyathon3 scraper/cazy_webscraper.py``.
