using System.IO;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.EventSystems;
using System.Drawing;
using System.Drawing.Imaging;
using SFB;

namespace FisheyeFriend {
    public class IO : MonoBehaviour {


        public Button openImageButton; // "Open Image" button

        public Button saveImageButton; // "Save" button

        public Button saveSharpImageButton; // "Save sharp (slow)" button

        public RawImage originalImage; // Left RawImage

        public RawImage processedImage; // Right RawImage

        public ScrollRect scrollRect; // ScrollRect for side-by-side images

        // Store the original loaded image
        private Texture2D originalTexture;

        // Store the resized image
        private Texture2D resizedTexture;

        private Vector2 screenSize;

        public GpuWarpUnity gpuWarp;


        void Start () {
            openImageButton.onClick.AddListener(OnOpenImageClicked);
            saveImageButton.onClick.AddListener(OnSaveImageClicked);
            saveImageButton.interactable = false; // Disable until an image is loaded
            saveSharpImageButton.onClick.AddListener(OnSaveSharpImageClicked);
            saveSharpImageButton.interactable = false; // Disable until an image is loaded
            screenSize = new Vector2(Screen.width, Screen.height);
            PositionButtons();

        }

        void PositionButtons () {
            Vector2 layoutCurr = new Vector2(-Screen.width/2, Screen.height/2);
            RectTransform currRect = openImageButton.GetComponent<RectTransform>();
            currRect.anchoredPosition = new Vector3(layoutCurr.x+ currRect.rect.width / 2, layoutCurr.y - currRect.rect.height / 2, openImageButton.transform.position.z);
            layoutCurr = new Vector2(currRect.anchoredPosition.x, Screen.height / 2);
            currRect = saveImageButton.GetComponent<RectTransform>();
            currRect.anchoredPosition = new Vector3(layoutCurr.x + currRect.rect.width, layoutCurr.y - currRect.rect.height / 2, saveImageButton.transform.position.z);
            layoutCurr = new Vector2(currRect.anchoredPosition.x, Screen.height / 2);
            currRect = saveSharpImageButton.GetComponent<RectTransform>();
            currRect.anchoredPosition = new Vector3(layoutCurr.x + currRect.rect.width, layoutCurr.y - currRect.rect.height / 2, saveSharpImageButton.transform.position.z);

        }

        void Update () {

            if (Screen.width != screenSize.x || Screen.height != screenSize.y) {
                PositionButtons();
            }
           if(Input.GetKeyDown(KeyCode.Tab)) {
                SwapImages();
            }
           if (Input.GetKeyDown(KeyCode.Escape)) {
                Application.Quit();
            }


           if(Input.GetKey(KeyCode.W)) {
                originalImage.transform.position -= new Vector3(0, 1* originalImage.transform.localScale.y, 0);
                processedImage.transform.position -= new Vector3(0, 1 * processedImage.transform.localScale.y, 0);
            }
            if (Input.GetKey(KeyCode.S)) {
                originalImage.transform.position += new Vector3(0, 1 * originalImage.transform.localScale.y, 0);
                processedImage.transform.position += new Vector3(0, 1 * processedImage.transform.localScale.y, 0);
            }

            if (Input.GetKey(KeyCode.A)) {
                originalImage.transform.position += new Vector3(1 * originalImage.transform.localScale.x, 0, 0);
                processedImage.transform.position += new Vector3(1 * processedImage.transform.localScale.x, 0, 0);
            }
            if (Input.GetKey(KeyCode.D)) {
                originalImage.transform.position -= new Vector3(1 * originalImage.transform.localScale.x, 0, 0);
                processedImage.transform.position -= new Vector3(1 * processedImage.transform.localScale.x, 0, 0);
            }

            if (Input.mouseScrollDelta.y != 0f) {
                float newScale = originalImage.transform.localScale.y+Input.mouseScrollDelta.y;
                if(newScale >0) {
                    originalImage.transform.localScale = new Vector3(newScale, newScale, 1);
                    processedImage.transform.localScale = new Vector3(newScale, newScale, 1);
                }
            }

            if (Input.GetKey(KeyCode.Minus)) {
                float newScale = originalImage.transform.localScale.y - originalImage.transform.localScale.y / 20f;
                if (newScale > 0) {
                    originalImage.transform.localScale = new Vector3(newScale, newScale, 1);
                    processedImage.transform.localScale = new Vector3(newScale, newScale, 1);
                }
            }

            if (Input.GetKey(KeyCode.Equals)) {
                float newScale = originalImage.transform.localScale.y + originalImage.transform.localScale.y / 20f;
                if (newScale > 0) {
                    originalImage.transform.localScale = new Vector3(newScale, newScale, 1);
                    processedImage.transform.localScale = new Vector3(newScale, newScale, 1);
                }
            }
        }

        void SwapImages () {
            if(!originalImage.gameObject.activeSelf && !processedImage.gameObject.activeSelf) {
                processedImage.gameObject.SetActive(true);
            } else if (originalImage.gameObject.activeSelf) {
                originalImage.gameObject.SetActive(false);
                processedImage.gameObject.SetActive(true);
            } else {
                originalImage.gameObject.SetActive(true);
                processedImage.gameObject.SetActive(false);
            }
        }

        void OnOpenImageClicked () {
            string path = ShowFilePicker();

            if (!string.IsNullOrEmpty(path) && File.Exists(path)) {
                //byte[] fileData = File.ReadAllBytes(path);
                //  originalTexture = new Texture2D(2, 2);

                Bitmap source = new Bitmap(path); // your loaded Bitmap
                Texture2D tex = BitmapToTexture2D(source);
                originalImage.texture = tex;

                //  if (originalTexture.LoadImage(fileData)) {

                //   byte[] sourceBytes = ((Texture2D)originalTexture.texture).EncodeToPNG();

                Texture2D wrarped = gpuWarp.WarpImage(source);
                   /* Bitmap warpedBM;
                    using (var warpMS = new MemoryStream(wrarped.EncodeToPNG())) {
                        warpedBM = new Bitmap(warpMS);
                    }
                    warpedBM.Save("D:\\tmp\\wallUnityGPULanc.png", ImageFormat.Png);
                   */
                    processedImage.texture = wrarped;

                // Enable the save button
                saveImageButton.interactable = true;
                saveSharpImageButton.interactable = true;
                SwapImages();
             //   }
            }
        }

        void OnSaveImageClicked () {
            string path = GetSavePath();
            if (!string.IsNullOrEmpty(path)) {
                TextureToBitmap((Texture2D)processedImage.texture).Save(path, ImageFormat.Png);
            }
        }

        string GetSavePath () {
            string path = ShowSavePicker();
            if (!string.IsNullOrEmpty(path)) {
                if (path.EndsWith(".png") == false) {
                    path += ".png";
                }
            }
            return path;
        }

        void OnSaveSharpImageClicked () {
            string path = GetSavePath();
            if (!string.IsNullOrEmpty(path)) {
                
                Bitmap sharpBM = new MapWarper().WarpImage(TextureToBitmap((Texture2D)originalImage.texture));
                sharpBM.Save(path, ImageFormat.Png);
            }
        }

        public Bitmap TextureToBitmap(Texture2D texture) {
            // Stashing an alternate way to convert
            /* byte[] bytes = texture.EncodeToPNG();
             using (var ms = new MemoryStream(bytes)) {
                 return new Bitmap(ms);
             }*/

            Bitmap bm;
            using (var warpMS = new MemoryStream(texture.EncodeToPNG())) {
                bm = new Bitmap(warpMS);
            }
            return bm;
        }


        public static Texture2D BitmapToTexture2D (Bitmap bitmap) {
            using (var memory = new MemoryStream()) {
                bitmap.Save(memory, System.Drawing.Imaging.ImageFormat.Png);
                memory.Position = 0;
                byte[] bytes = memory.ToArray();
                Texture2D tex = new Texture2D(2, 2);
                tex.LoadImage(bytes);
                return tex;
            }
        }
        /**
         resizedTexture = ResizeTextureMaintainingAspect(originalTexture, 1920/2, 1080);

                    // Display the resized image on both sides
                    leftImage.texture = resizedTexture;
                    rightImage.texture = Instantiate(resizedTexture);

                    // Adjust the size of the RawImages to fit the resized texture
                    Vector2 size = new Vector2(resizedTexture.width, resizedTexture.height);
                    leftImage.rectTransform.sizeDelta = size;
                    rightImage.rectTransform.sizeDelta = size;

         Bitmap rightBM;
                    using (var ms = new MemoryStream(pngData)) {
                        rightBM = new Bitmap(ms);
                    }
        */


        // Resize the texture maintaining aspect ratio
        Texture2D ResizeTextureMaintainingAspect (Texture2D source, int maxWidth, int maxHeight) {
            float ratio = Mathf.Min((float)maxWidth / source.width, (float)maxHeight / source.height, 1f);
            int newWidth = Mathf.RoundToInt(source.width * ratio);
            int newHeight = Mathf.RoundToInt(source.height * ratio);

            Texture2D result = new Texture2D(newWidth, newHeight, source.format, false);

            for (int y = 0; y < newHeight; y++) {
                for (int x = 0; x < newWidth; x++) {
                    float u = x / (float)(newWidth - 1);
                    float v = y / (float)(newHeight - 1);
                    UnityEngine.Color color = source.GetPixelBilinear(u, v);
                    result.SetPixel(x, y, color);
                }
            }
            result.Apply();
            return result;
        }


        string ShowFilePicker () {
            var paths = StandaloneFileBrowser.OpenFilePanel("Open File", "", "", false);
            return paths[0];
        }

        string ShowSavePicker () {
            var path = StandaloneFileBrowser.SaveFilePanel("Save File", "", "", "");
            return path;
        }




    }
}