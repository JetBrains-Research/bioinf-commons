[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
[![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Biolabs_BioinfCommons)/statusIcon.svg)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Biolabs_BioinfCommons&guest=1)

Bioinf-commons
==============
Bioinformatics library in [Kotlin](https://kotlinlang.org) by [JetBrains Research](https://research.jetbrains.org/groups/biolabs). 

* APIs for working with Genome, Sequence, Genes, Ontologies etc.
* BAM, BED, FAST* formats support
* Pandas-like data frames
* Statistics utilities including distributions mixtures, hmms and hypothesis testing
* Utilities for computations caching, logging, parallel execution with cancelling, etc.

Installation
------------

The latest version of `bioinf-commons` is available on [jCenter](https://bintray.com/bintray/jcenter). If you're using
Gradle just add the following to your `build.gradle`:

```gradle
repositories {
    jcenter()
}

dependencies {
    compile 'org.jetbrains.bio:bioinf-commons:0.0.8'
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