[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
[![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Biolabs_BioinfCommons)/statusIcon.svg)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Biolabs_BioinfCommons&guest=1)

Bioinf-commons
==============
Bioinformatics library in [Kotlin](https://kotlinlang.org).

Contents
--------

* `org.jetbrains.bio.dataframe` - Pandas like dataframe
* `org.jetbrains.bio.experiment` - Named computation - experiment and resources configuration
* `org.jetbrains.bio.genome` - APIs for working with Genome, Sequence, Genes, Ontologies etc.
    * `containers` - Genome location sets API
    * `coverage` - Genome coverage API, paired/single end, fragment size estimation
    * `data` - Describe any kind of dataset resources including replicates
    * `format` - Bam (including Bisulfite sequencing), Bed, Fasta, Fastq, 2bit formats support
    * `methylome` - API to work with methylomes - filtration, statistics, aggregations etc.
    * `query` - Named functions - queries with caching capabilities
    * `sampling` - API for genomic sampling - sequencies, locations, etc.
    * `sequence` - Genome sequence API<br/>
    Also: Biomart, Ensembl, UCSC support, Genomes and Genes annotations
* `org.jetbrains.bio.statistics` - Statistics utilities including distributions mixtures, hmms and hypothesis testing
* `org.jetbrains.bio.util` - Cancellable computations, progress reporters, logging utilities, and other utils

Installation
------------

The latest version of `bioinf-commons` is available on [jCenter](https://bintray.com/bintray/jcenter). If you're using
Gradle just add the following to your `build.gradle`:

```gradle
repositories {
    jcenter()
}

dependencies {
    compile 'org.jetbrains.bio:bioinf-commons:0.0.9'
}
```

Publishing
----------

You can publish a new release with a one-liner

```bash
./gradlew clean assemble test generatePomFileForMavenJavaPublication bintrayUpload
```

Make sure to set Bintray credentials (see API key section
[here](https://bintray.com/profile/edit)) in `$HOME/.gradle/gradle.properties`.

```
$ cat $HOME/.gradle/gradle.properties
bintrayUser=CHANGEME
bintrayKey=CHANGEME
```
