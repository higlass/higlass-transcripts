import { expect } from "chai";
import register from "higlass-register";

import { getTrackObjectFromHGC } from "higlass";

import {
  mountHGComponent,
  removeHGComponent,
} from "./utils/test-helpers";

import viewConf from "./view-configs/simple-track-aa";

import TranscriptsTrack from "../src/scripts/TranscriptsTrack";

register({
  name: "TranscriptsTrack",
  track: TranscriptsTrack,
  config: TranscriptsTrack.config,
});



describe("Transcripts tests", () => {
  let hgc = null;
  let div = null;

  beforeAll((done) => {
    [div, hgc] = mountHGComponent(div, hgc, viewConf, done);
  });

  it("tests that basic infos are correct", (done) => {

    const trackObj = getTrackObjectFromHGC(
      hgc.instance(),
      viewConf.views[0].uid,
      viewConf.views[0].tracks.top[0].uid
    );

    //const tile = trackObj.visibleAndFetchedTiles()[2];

    expect(trackObj.trackHeight).to.equal(106);
    expect(Object.keys(trackObj.transcriptInfo).length).to.equal(5);

    const transcriptInfo = trackObj.transcriptInfo["ENST00000270722.10_chr1_3069203_3438621"];
    expect(transcriptInfo.transcriptName).to.equal("PRDM16-201");
    expect(transcriptInfo.txStart).to.equal(3069202);
    expect(transcriptInfo.displayOrder).to.equal(3);

    // Check position info
    const posInfo = trackObj.transcriptPositionInfo[3][0]
    expect(posInfo[1]).to.equal(3438621);
    expect(trackObj.areTranscriptsHidden).to.equal(false);
    done();

  });

  afterAll(() => {
    removeHGComponent(div);
  });
});


