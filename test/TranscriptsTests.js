import { expect } from "chai";
import register from "higlass-register";

import FetchMockHelper from "./utils/FetchMockHelper";

import { HiGlassComponent, getTrackObjectFromHGC } from "higlass";

import {
  waitForDataLoaded,
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

describe("Transcripts", () => {
  const fetchMockHelper = new FetchMockHelper("", "TranscriptTests");

  beforeAll(async () => {
    await fetchMockHelper.activateFetchMock();
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

      const transcriptInfo = trackObj.transcriptInfo["ENST00000270722.9_chr1_3069211_3438621"];
      expect(transcriptInfo.transcriptName).to.equal("PRDM16-201");
      expect(transcriptInfo.txStart).to.equal(3069210);
      expect(transcriptInfo.displayOrder).to.equal(4);

      // Check position info
      const posInfo = trackObj.transcriptPositionInfo[3][0]
      expect(posInfo[1]).to.equal(3434342);

      expect(trackObj.areTranscriptsHidden).to.equal(false);

      // Check vertical positioning of exon
      
      // const exonPos1 = tile.allExonsForSVG[4].rect[1];
      // const exonPos2 = tile.allExonsForSVG[4].rect[5];
      // expect(exonPos1).to.equal(41);
      // expect(exonPos2).to.equal(47);

      done();

      // setTimeout(() => {

      // }, 1000);


      
    });

    afterAll(() => {
      removeHGComponent(div);
    });
  });

  afterAll(async () => {
    await fetchMockHelper.storeDataAndResetFetchMock();
  });
});
