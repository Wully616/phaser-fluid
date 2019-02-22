import FluidGrid from '../objects/liquid/FluidGrid.js'
import game from '../main/game.js'

export default class SceneMain extends Phaser.Scene {
    constructor() {
        super({ key: "SceneMain" });
        this.blocksize = 8;
        
    }

    preload(){

    }

    create(){
        this.text = this.add.text(10, 10, '', { fill: '#ff0000' }).setDepth(1);
        this.fluid = new FluidGrid(this,0,0,"fluid",game.config.width,game.config.height,this.blocksize,0.0001,0.0000001,0.01);
        this.input.mouse.disableContextMenu();

        
        this.input.on('pointerdown', function (pointer) {
            if (pointer.rightButtonDown()){
                this.rightClick = {
                    downX: pointer.x,
                    downY: pointer.y
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
        //this.fluid.addDensity(game.config.width/this.blocksize*0.5,game.config.height/this.blocksize*0.5,100);
        //this.fluid.addDensity(game.config.width/this.blocksize*0.5+1,game.config.height/this.blocksize*0.5,100);
        //this.fluid.addDensity(game.config.width/this.blocksize*0.5-1,game.config.height/this.blocksize*0.5,100);
        //this.fluid.addDensity(game.config.width/this.blocksize*0.5,game.config.height/this.blocksize*0.5-1,100);
        //this.fluid.addDensity(game.config.width/this.blocksize*0.5,game.config.height/this.blocksize*0.5+1,100);
        //this.fluid.addVelocity(game.config.width/this.blocksize*0.5,game.config.height/this.blocksize*0.5,0,-20); 
        
        /*
        this.text.setText([
            'x: ' + pointer.worldX,
            'y: ' + pointer.worldY,
            'gridX: ' + gridX,
            'gridY: ' + gridY,
            'gridDensity: ' + this.fluid.getDensity(gridX, gridY),
            'gridVelocityX: ' + this.fluid.getVelocityX(gridX, gridY),
            'gridVelocityY: ' + this.fluid.getVelocityY(gridX, gridY),
            'isDown: ' + pointer.isDown,
            'rightButtonDown: ' + pointer.rightButtonDown()
        ]);
        */   
        

        if(pointer.leftButtonDown()){          
            var gridX = Math.round(pointer.x / this.blocksize);
            var gridY = Math.round(pointer.y / this.blocksize);  
            this.fluid.addDensity(gridX,gridY,1000); 
            this.fluid.addDensity(gridX+1,gridY,1000); 
            this.fluid.addDensity(gridX-1,gridY,1000); 
            this.fluid.addDensity(gridX,gridY+1,1000); 
            this.fluid.addDensity(gridX,gridY-1,1000); 

            this.fluid.addDensity(gridX+1,gridY+1,1000); 
            this.fluid.addDensity(gridX-1,gridY-1,1000); 
            this.fluid.addDensity(gridX-1,gridY+1,1000); 
            this.fluid.addDensity(gridX+1,gridY-1,1000);

        }
        if(pointer.rightButtonDown()){
            var gridX = Math.round(this.rightClick.downX / this.blocksize);
            var gridY = Math.round(this.rightClick.downY / this.blocksize);

            var vX = pointer.x - this.rightClick.downX;
            var vY = pointer.y - this.rightClick.downY;
            this.fluid.addVelocity(gridX,gridY,vX,vY);  
        }
    }
}