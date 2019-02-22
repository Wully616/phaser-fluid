import SceneMainMenu from '../scene/SceneMainMenu.js'
import SceneMain from '../scene/SceneMain.js'

var config = {
    type: Phaser.WEBGL,
    width: 512,
    height: 512,
    
    physics: {
        default: "arcade",
        arcade: {
            gravity: { x: 0, y: 0 }
      }
    },
    scene: [
        SceneMainMenu,
        SceneMain
    ],
    pixelArt: true,
    roundPixels: true
};


var game = new Phaser.Game(config);
export default game;