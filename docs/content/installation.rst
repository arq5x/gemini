############
Installation
############

Automated installation
======================

GEMINI contains an automated installation script which installs
GEMINI along with required Python dependencies, third party software
and data files. 

.. code-block:: bash

    $ wget https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py

If ``wget`` isn't available, you can use ``curl`` instead:

.. code-block:: bash

    $ curl \
    https://raw.github.com/arq5x/gemini/master/gemini/scripts/gemini_install.py \
    > gemini_install.py

Once you have downloaded the above installation script, you can proceed as follows:

.. code-block:: bash

    $ python gemini_install.py /usr/local /usr/local/share/gemini
    $ export PATH=$PATH:/usr/local/gemini/bin
    # it would be wise to add the above line to your ``.bashrc`` or ``.bash_profile``

This installs the GEMINI executable as ``/usr/local/bin/gemini``,
other required third party dependencies in ``/usr/local/bin``, and
associated data files in ``/usr/local/share/gemini``.  **Please note that this
is merely an example: one can place the GEMINI executables and annotation files
in any directories one wishes.**

.. note::

  The automated installtion script typically requires ~30 minutes, primarily
  owing to the time required to download the GEMINI genome annotation files.
  Also, please note that the annotation files **requires ~15Gb of storage**, 
  so please ensure that the data directory (``/usr/local/share/gemini`` in 
  the example above) has sufficient space.

.. tip::

  **Some tips and tricks for installation issues:**

  1. Some older versions of wget have certificate problems with GitHub
     files. If you run into this problem, you can alternatively download
     the install script using``wget --no-check-certificates`` or ``curl -O``.

  2. The installation script is idempotent and you can re-run it multiple
     times without any issues. If you experience internet connectivity or
     other transient errors during installation, a re-run can often solve
     the problem (fingers crossed).

  3. If you're installing behind a proxy you'll need to set proxy information in
     a `condarc file`_ and also set ``all_proxy``, ``http_proxy`` and ``https_proxy``
     in your ``.bashrc`` file. See `this mailing list discussion`_ for more information.

.. _condarc file: http://docs.continuum.io/conda/intro.html#configuration
.. _this mailing list discussion: https://groups.google.com/d/msg/gemini-variation/WMQiMDeW6aE/6oa8nS1NorUJ


Dependencies
-------------------------------
The installer requires:
  
  - Python 2.7.x
  - git
  - a working C / C++ compiler such as gcc
  - zlib (including headers)

These dependencies already exist on most UNIX/LINUX/OS X machines.  However,
on minimalist machines **such as fresh Amazon Cloud instances**, one may have
to install these dependencies before running the automated installer. In the
case of Amazon instances, the following command should take care of all of the
above dependencies:

.. code-block:: bash

  $ sudo yum -y install python27 git gcc gcc-c++ zlib-devel


Installing without root access.
-------------------------------
As many users do not have root or sudo access, the automated installer also 
has options to install in "non-root" environments:

.. code-block:: bash

  $ python gemini_install.py ~/gemini ~/gemini --nosudo


Updating your PATH to find the GEMINI executable
-------------------------------------------------
At this point, you will have a self-contained installation of GEMINI, 
including both the software and its associated genome annotations. However,
if you have done a custom install in a "non-root" enviornment, you will
first need to update your ``PATH`` environment variable to include the path
to the bin directory that you just created by running the automated installer.

For example, if, as above, you placed you custom install in ``~/gemini``, you
would need to update your ``PATH`` as follows. It would be wise to also add this
to your ``.bashrc`` or ``.bash_profile``:

.. code-block:: bash

    $ export PATH=$PATH:~/gemini/bin

Note that this change will only last for the life of your current terminal 
session.  To make this more permanent, update your ``.bash_profile`` so that
this change is made each time you login.

If successful, you should be able to run the following command from anywhere
on your system:

.. code-block:: bash

   $ gemini -v
   gemini 0.3.0b


Running unit tests
-------------------------------------------------
If successfully installed, you should be able to change
directories into the ``gemini`` directory within the 
directory into which you installed the GEMINI source code
and run a script of unit tests that will ensure that GEMINI
is running appropriately on your system. For example, if, as 
above, you installed the GEMINI executables to ``~/gemini``, you
would issue the following commands to run the unit tests:

.. code-block:: bash

   $ cd ~/gemini/gemini
   $ bash master-test.sh 

Updating the GEMINI executables and annotations
-------------------------------------------------
Once installed with the automated installer, it is easy to upgrade the
GEMINI programs and associated data files to the latest released 
version with:

.. code-block:: bash

    $ gemini update

There are also flags to update to the latest development version of the code or
to only update data files without updating the source:

.. code-block:: bash

    $ gemini update --devel
    $ gemini update --dataonly


Software dependencies
=====================
GEMINI depends upon several widely-used genomics command line software as well
as multiple Python packages.  We recognize that the dependency stack is quite
deep and are working on ways to minimize dependencies in the interest of the
most streamlined installation process possible.  Nonetheless, the following are
core dependencies:

    1. Python 2.7.x
    2. `grabix <https://github.com/arq5x/grabix>`_
    3. `samtools <http://sourceforge.net/projects/samtools/files/>`_
    4. `tabix <http://sourceforge.net/projects/samtools/files/>`_
    5. `bedtools <https://code.google.com/p/bedtools/>`_
    6. `pybedtools <http://pythonhosted.org/pybedtools/main.html#installing-pybedtools>`_

Manual installation
=====================
Once the above dependencies have been installed, one can begin installing
``GEMINI`` itself. To install you should download the latest source code from
GitHub, either by going to::

    http://github.com/arq5x/gemini

and clicking on "Downloads", or by cloning the git repository with:

.. code-block:: bash

    $ git clone https://github.com/arq5x/gemini.git

Once you have the source code, run:

.. code-block:: bash

    $ cd gemini
    $ sudo python setup.py install

to install it. If you don't have permission to install it in the default
directory, you can simply build the source in-place and use the package
from the git repository:

.. code-block:: bash

    $ python setup.py build_ext --inplace


Installing annotation files
===========================
One of the more appealing features in ``GEMINI`` is that it automatically
annotates variants in a VCF file with several genome annotations.  However,
you must first install these data files on your system. It's easy enough ---
you just need to run the following script and tell it in which what full path
you'd like to install the necessary data files. The recommended path is
``/usr/local/share``, but you can install the data files wherever you want.

.. code-block:: bash

    $ python gemini/install-data.py /usr/local/share/
    
Using previously installed annotation files
===============================================================
If you have installed GEMINI with the annotation files on a server and you can NFS mount
the annotation files, you can tell a local install of GEMINI where those annotation files
are by making the file ~/.gemini/gemini-config.yaml::

	annotation_dir: /path/to/nfs_mounted/gemini/data
	versions:
  	  GRCh37-gms-mappability.vcf.gz: 2
  	  hg19.rmsk.bed.gz: 2

Running the testing suite
===========================
GEMINI comes with a full test suite to make sure that everything has installed
correctly on your system.  We **strongly** encourage you to run these tests.

.. code-block:: bash

    $ bash master-test.sh


---------------------------
Functional annotation tools
---------------------------
`GEMINI` depends upon external tools to predict the functional consequence of variants in a VCF file.
We currently support annotations produced by both `SnpEff <http://snpeff.sourceforge.net/>`_
and `VEP <http://www.ensembl.org/info/docs/variation/vep/index.html>`_.
Recommended instructions for annotating existing VCF files with these tools are available here.
In addition, we have attempted to standardize the terms used to describe the functional consequence of a given variant,
as each annotation tool uses different vocabulary.


The variant consequence columns in the variant table are populated either by `snpEff` or `VEP` as defined by the user using the `-t` option while running gemini load
(To populate these columns the input VCF file should have been annotated either by `snpEff` or `VEP`):

.. code-block:: bash

	$ gemini load -v my.vcf -t VEP -d my.db
	$ gemini load -v my.vcf -t snpEFF -d my.db

By default the following columns in the variant table would be set to null:

* anno_id
* gene
* affected_gene
* affected_transcript
* affected_exon
* is_exonic
* is_lof
* is_coding
* codon_change
* aa_change
* aa_length
* biotype
* most_severe_impact
* impact_severity
* polyphen_pred
* polyphen_score
* sift_pred
* sift_score

Impacts
.......
The table below shows the alternate `GEMINI` terms for the consequences from `snpEff` and `VEP`, for SQL queries.
The last column represents the severity terms associated with the impacts:

=============================================       ===================================    =====================================================     ================
Gemini terms                                        snpEff terms                           VEP terms                                                 Impact severity
=============================================       ===================================    =====================================================     ================
splice_acceptor                                     SPLICE_SITE_ACCEPTOR                   splice_acceptor_variant                                   HIGH
splice_donor                                        SPLICE_SITE_DONOR                      splice_donor_variant                                      HIGH
stop_gain                                           STOP_GAINED                            stop_gained                                               HIGH
stop_loss                                           STOP_LOST                              stop_lost                                                 HIGH
frame_shift                                         FRAME_SHIFT                            frameshift_variant                                        HIGH
start_loss                                          START_LOST                             null                                                      HIGH
exon_deleted                                        EXON_DELETED                           null                                                      HIGH
non_synonymous_start                                NON_SYNONYMOUS_START                   null                                                      HIGH
transcript_codon_change                             null                                   initiator_codon_variant                                   HIGH
non_syn_coding                                      NON_SYNONYMOUS_CODING                  missense_variant                                          MED
inframe_codon_gain                                  CODON_INSERTION                        inframe_insertion                                         MED
inframe_codon_loss                                  CODON_DELETION                         inframe_deletion                                          MED
inframe_codon_change                                CODON_CHANGE                           null                                                      MED
codon_change_del                                    CODON_CHANGE_PLUS_CODON_DELETION       null                                                      MED
codon_change_ins                                    CODON_CHANGE_PLUS_CODON_INSERTION      null                                                      MED
UTR_5_del                                           UTR_5_DELETED                          null                                                      MED
UTR_3_del                                           UTR_3_DELETED                          null                                                      MED
other_splice_variant                                null                                   splice_region_variant                                     MED
mature_miRNA                                        null                                   mature_miRNA_variant                                      MED
regulatory_region                                   null                                   regulatory_region_variant                                 MED
TF_binding_site                                     null                                   TF_binding_site_variant                                   MED
regulatory_region_ablation                          null                                   regulatory_region_ablation                                MED
regulatory_region_amplification                     null                                   regulatory_region_amplification                           MED
TFBS_ablation                                       null                                   TFBS_ablation                                             MED
TFBS_amplification                                  null                                   TFBS_amplification                                        MED
synonymous_stop                                     SYNONYMOUS_STOP                        stop_retained_variant                                     LOW
synonymous_coding                                   SYNONYMOUS_CODING                      synonymous_variant                                        LOW
UTR_5_prime                                         UTR_5_PRIME                            5_prime_UTR_variant                                       LOW
UTR_3_prime                                         UTR_3_PRIME                            3_prime_UTR_variant                                       LOW
intron                                              INTRON                                 intron_variant                                            LOW
CDS                                                 CDS                                    coding_sequence_variant                                   LOW
upstream                                            UPSTREAM                               upstream_gene_variant					                 LOW
downstream                                          DOWNSTREAM                             downstream_gene_variant                                   LOW
intergenic                                          INTERGENIC, INTERGENIC_CONSERVED       intergenic_variant                                        LOW
intragenic                                          INTRAGENIC                             null                                                      LOW
gene                                                GENE                                   null                                                      LOW
transcript                                          TRANSCRIPT                             null                                                      LOW
exon                                                EXON                                   null                                                      LOW
start_gain                                          START_GAINED                           null                                                      LOW
synonymous_start                                    SYNONYMOUS_START                       null                                                      LOW
intron_conserved                                    INTRON_CONSERVED                       null                                                      LOW
nc_transcript                                       null                                   nc_transcript_variant                                     LOW
NMD_transcript                                      null                                   NMD_transcript_variant                                    LOW
incomplete_terminal_codon                           null                                   incomplete_terminal_codon_variant                         LOW
nc_exon                                             null                                   non_coding_exon_variant                                   LOW
transcript_ablation                                 null                                   transcript_ablation                                       LOW
transcript_amplification                            null                                   transcript_amplification                                  LOW
feature elongation                                  null                                   feature elongation                                        LOW
feature truncation                                  null                                   feature truncation                                        LOW
=============================================       ===================================    =====================================================     ================
*Note: "null" refers to the absence of the corresponding term in the alternate database*
