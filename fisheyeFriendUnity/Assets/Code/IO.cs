using System.IO;
using UnityEngine;
using UnityEngine.UI;
using UnityEngine.EventSystems;
using System.Drawing;
using System.Drawing.Imaging;
using SFB;

namespace FisheyeFriend {
    public class IO : MonoBehaviour {


        // Assign in the Inspector: Button for opening file dialog
        public Button openImageButton; // "Open Image" button

        // Assign in the Inspector: Button for saving the original image
        public Button saveImageButton; // "Save" button

        // Assign in the Inspector: RawImage for displaying the resized image (left)
        public RawImage originalImage; // Left RawImage

        // Assign in the Inspector: RawImage for displaying the resized image (right)
        public RawImage processedImage; // Right RawImage

        // Assign in the Inspector: ScrollRect containing both images
        public ScrollRect scrollRect; // ScrollRect for side-by-side images

        // Store the original loaded image
        private Texture2D originalTexture;

        // Store the resized image
        private Texture2D resizedTexture;



        public GpuWarpUnity gpuWarp;


        void Start () {
            openImageButton.onClick.AddListener(OnOpenImageClicked);
            saveImageButton.onClick.AddListener(OnSaveImageClicked);
            saveImageButton.interactable = false; // Disable until an image is loaded

            Vector3 leftEdge = Camera.main.ScreenToWorldPoint(new Vector3(0, Screen.height, openImageButton.transform.position.z));
            Vector3 topEdge = Camera.main.ScreenToWorldPoint(new Vector3(Screen.width, 0, openImageButton.transform.position.z));

           // openImageButton.transform.position = new Vector3(leftEdge.x + openImageButton.GetComponent<RectTransform>().rect.width / 2, topEdge.y - openImageButton.GetComponent<RectTransform>().rect.height / 2, openImageButton.transform.position.z);
        }

        void Update () {
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
                originalImage.transform.localScale += new Vector3(Input.mouseScrollDelta.y, Input.mouseScrollDelta.y, 0);
                processedImage.transform.localScale += new Vector3(Input.mouseScrollDelta.y, Input.mouseScrollDelta.y, 0);
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
            Debug.Log("Open Image button clicked.");
            // Use a file picker plugin or a custom file dialog for Unity.
            // For demonstration, this is a placeholder path.
            string path = ShowFilePicker(); // Implement this with a plugin or native dialog

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
                SwapImages();
             //   }
            }
        }

        void OnSaveImageClicked () {
            TextureToBitmap((Texture2D)processedImage.texture).Save("D:\\tmp\\wallUnityGPULancSave.png", ImageFormat.Png);
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





    }
}