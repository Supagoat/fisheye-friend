using System.Collections;
using System.Collections.Generic;
using UnityEngine;

namespace FisheyeFriend {
    public class UnityLink : MonoBehaviour {
        private int ct = 0;
        void Start () {

        }

        // Update is called once per frame
        void Update () {
            ct++;
            if (ct == 1) {

                WarpConfig config = new WarpConfig();
                config.mapPath = "D:\\tmp\\hemi_pairs_596.json";
                config.srcPath = "D:\\tmp\\wallFish.jpg";
                config.outPath = "D:\\tmp\\wallUnityOut.png";

                MapWarper warper = new MapWarper();
                warper.WarpImage(config);
            }
        }
    }
}
