# HiGlass Transcripts Track - Under development

Disply gene transcripts in HiGlass!

Zoomed out:

![Transcripts track](https://aveit.s3.amazonaws.com/higlass/static/higlass-transcripts-zoomed-out.png)

Zoomed in:

![Transcripts track](https://aveit.s3.amazonaws.com/higlass/static/higlass-transcripts-zoomed-in.png)

**Note**: This is the source code for the transcripts track only! You might want to check out the following repositories as well:

- HiGlass viewer: https://github.com/higlass/higlass
- HiGlass server: https://github.com/higlass/higlass-server
- HiGlass docker: https://github.com/higlass/higlass-docker

## Installation
 
```
npm install higlass-transcripts
```

## Data preparation

To extract trancript data from a Gencode GTF file, the following script can be used (make sure to adjust the file names in the script)
```
python /scripts/extract_transcript_data.py
```

To create an aggregated `beddb` file from that data, you can use the script
```
python /scripts/aggregate_transcripts.py
```

To ingest the data into higlass-server:
```
python manage.py ingest_tileset --filename data/transcripts.beddb --filetype beddb --datatype gene-annotation --uid aweseome_transcripts
```


## Usage

The live script can be found at:

- https://unpkg.com/higlass-transcripts/dist/higlass-transcripts.js

### Client

1. Make sure you load this track prior to `hglib.js`. For example:

```
<script src="/higlass-transcripts.js"></script>
<script src="hglib.js"></script>
<script>
  ...
</script>
```

### Options
The following options are available:
```
{
  "server": "http://localhost:8001/api/v1",
  "tilesetUid": "awesome_transcripts",
  "uid": "awesome_transcripts_uid",
  "type": "horizontal-transcripts",
  "options": {
    "fontSize": 9, // font size for labels and amino acids (if available)
    "fontFamily": "Helvetica",
    "labelFontColor": "#222222",
    "labelBackgroundColor": "#e9e9e9",
    "plusStrandColor": "#bdbfff", // color of coding parts of the exon on the plus strand
    "minusStrandColor": "#fabec2", // color of coding parts of the exon on the negative strand
    "utrColor": "#C0EAAF", // color of untranslated regions of the exons
    "transcriptHeight": 12, // height of the transcripts
    "transcriptSpacing": 2, // space in between the transcripts
    "name": "Gene transcripts",
    "maxTexts": 50, // Maximum number of labels shown on the screen
    "showToggleTranscriptsButton": true, // If the "Show fewer transcripts"/"Show more transcripts" is shown
    "trackHeightAdjustment": "automatic", // if "automatic", the height of the track is adjusted to the number of visible transcripts.
    "startCollapsed": false, // if true, only one transcript is shown
    "sequenceData": { // If this is set, transcribed amino acids are displayed when sufficiently zoomed in
      "type": "fasta",
      "fastaUrl": "https://aveit.s3.amazonaws.com/higlass/data/sequence/hg38.fa",
      "faiUrl": "https://aveit.s3.amazonaws.com/higlass/data/sequence/hg38.fa.fai",
      "chromSizesUrl": "https://aveit.s3.amazonaws.com/higlass/data/sequence/hg38.mod.chrom.sizes"
    },
  },
  "width": 768,
  "height": 40
}
```

### ECMAScript Modules (ESM)

We also build out ES modules for usage by applications who may need to import or use `higlass-transcripts` as a component.

Whenever there is a statement such as the following, assuming `higlass-transcripts` is in your node_modules folder:
```javascript
import { TranscriptsTrack } from 'higlass-transcripts';
```

Then TranscriptsTrack would automatically be imported from the `./es` directory (set via package.json's `"module"` value). 

## Support

For questions, please either open an issue or ask on the HiGlass Slack channel at http://bit.ly/higlass-slack

## Development

### Testing

To run the test suite:

```
npm run test-watch
```


### Installation

```bash
$ git clone https://github.com/higlass/higlass-transcripts.git
$ cd higlass-transcripts
$ npm install
```
If you have a local copy of higlass, you can then run this command in the higlass-transcripts directory:

```bash
npm link higlass
```

### Commands

 - **Developmental server**: `npm start`
 - **Production build**: `npm run build`
