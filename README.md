# Legacy STRAIGHT

The legacy-STRAIGHT is a collection of speech analysis, modification and resynthesis tools. 

## Installation

Set MATLAB path to "src" directory.

## Quick start

Paste the following code to MATLAB command window. It generates copy-synthesized output to the variable "syntheszed_signal".

    [x, fs] = audioread('vaiueo2d.wav'); 
    f0raw = MulticueF0v14(x,fs); 
    ap = exstraightAPind(x,fs,f0raw);
    n3sgram=exstraightspec(x,f0raw,fs);
    syntheszed_signal = exstraightsynth(f0raw,n3sgram,ap,fs);

For running this using GNU Octave, please load signal package.

```
pkg load signal
```

## Release note

* [July 19, 2018; Prerelease]
The "Quick start" example also runs properly on GNU Octave 4.4.0 on macOS High Sierra (10.13.6)

* [July 17, 2018: Prerelease]
Added documemts. The first release will be on July 24, 2018.

* [July 16, 2018: Prerelease]
This release is a copy of the latest version which was distributed by the first author (Hideki Kawahara) to academic communities. The version is named STRAIGHTV40_007. The last update was July 17, 2016. Kansai TLO has also licensed the legacy-STRAIGHT for commercial use. The licensees of the legacy-STRAIGHT agreed to make the legacy-STRAIGHT open to the public after July 15, 2018. 

## Acknowledgment

The legacy-STRAIGHT was supported by many coauthors, contributors, and funding agencies. 

***

Hideki Kawahara,
July 16, 2018 (start date)