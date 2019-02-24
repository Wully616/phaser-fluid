import FluidGrid from '../objects/liquid/FluidGrid.js'
import game from '../main/game.js'

export default class SceneMain extends Phaser.Scene {
    constructor() {
        super({ key: "SceneMain" });
        this.blocksize = 16;
        this.itemTypes = [
            "red","green","blue","addwall", "removewall"
        ];
        this.selectedItem = 0;
    }

    preload(){
        
    }

    create(){
        this.text = this.add.text(10, 10, '', { fill: '#ff0000' }).setDepth(1);
        this.fluid = new FluidGrid(this,0,0,"fluid",game.config.width,game.config.height,this.blocksize,0.000,0.0,0.1);
        this.input.mouse.disableContextMenu();

        
        
        this.input.on('pointerdown', function (pointer) {
            //gets the position the mouse was at when clicked
            if (pointer.rightButtonDown()){
                this.rightClick = {
                    downX: pointer.x,
                    downY: pointer.y
                }
            }
            if(pointer.middleButtonDown()){
                var length = this.itemTypes.length - 1;
                this.selectedItem++;
                if( this.selectedItem > length){
                    this.selectedItem = 0;
                }
            }
        }, this);

        

 
    }

    update(){
        this.fluid.fluidStep();
        this.fluid.render();
  
        var pointer = this.input.activePointer;
        var gridX = Math.floor(pointer.x / this.blocksize);
        var gridY = Math.floor(pointer.y / this.blocksize);

        
        this.text.setText([
            'x: ' + pointer.worldX,
            'y: ' + pointer.worldY,
            'r,g,b Density: ' + this.fluid.getDensity(gridX, gridY),
            
            'gridVelocityX: ' + this.fluid.getVelocityX(gridX, gridY),
            'gridVelocityY: ' + this.fluid.getVelocityY(gridX, gridY),
            'selected Item: ' + this.itemTypes[this.selectedItem]
            
        ]);
          
        
        
        if(pointer.leftButtonDown()){          
            var gridX = Math.round(pointer.x / this.blocksize);
            var gridY = Math.round(pointer.y / this.blocksize);  

            this.fluid.addItem(gridX,gridY,this.itemTypes[this.selectedItem], 1000);

        }
        if(pointer.rightButtonDown()){
            var gridX = Math.round(this.rightClick.downX / this.blocksize);
            var gridY = Math.round(this.rightClick.downY / this.blocksize);

            var vX = pointer.x - this.rightClick.downX;
            var vY = pointer.y - this.rightClick.downY;
            this.fluid.addVelocity(gridX,gridY,vX*0.1,vY*0.1);  
        }
    }
}