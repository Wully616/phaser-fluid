//https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
//http://www.cs.northwestern.edu/~sco590/cs140/assignment3.html
import Vector from '../../util/Vector2D.js'
export default class FluidGrid extends Phaser.GameObjects.GameObject {

    constructor(scene, x, y, key, gridWidth,gridHeight,blockSize, diffusion, viscosity, dt) {
        super(scene, x, y, key);
        this.scene = scene;
        this.scene.add.existing(this);

        this.gridWidth = gridWidth;
        this.gridHeight = gridHeight;
        this.blockSize = blockSize;
        //a debug grid
        this.createDebugGrid(this.scene,this.gridWidth,this.gridHeight,this.blockSize);
        
        this.iter = 1;
        //size of the grid
        this.size = gridWidth/blockSize;
        this.gridsize = this.size*this.size;
        this.N = this.size;
        //the timestamp, diffusion and viscosity of this fluidgrid
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        this.sR = new Array(this.gridsize).fill(0);
        this.sG = new Array(this.gridsize).fill(0);
        this.sB = new Array(this.gridsize).fill(0);
        this.densityR = new Array(this.gridsize).fill(0);
        this.densityG = new Array(this.gridsize).fill(0);
        this.densityB = new Array(this.gridsize).fill(0);
        this.Vx = new Array(this.gridsize).fill(0);
        this.Vx0 = new Array(this.gridsize).fill(0);
        this.Vy = new Array(this.gridsize).fill(0);
        this.Vy0 = new Array(this.gridsize).fill(0);

        this.solids = new Array(this.gridsize).fill(false);
        this.tmpar = []; //three temp arrays
        for(var ar = 0 ; ar < 3; ar ++){
            this.tmpar[ar] = new Array(this.gridsize).fill(0);
        }
        this.graphics = this.scene.add.graphics();

        this.wrap = true; // wrap grid bounds
        this.ixList = {};

        for (var i = -10; i < this.N +10 ; i++) {
            this.ixList[i] = {};
            for (var j = -10; j < this.N+10; j++) { 
                
                this.ixList[i][j] = {
                    up: this.IXWrap(i,j - 1),
                    down: this.IXWrap(i,j + 1),
                    left: this.IXWrap(i-1 ,j),
                    right: this.IXWrap(i+1,j),
                    ix: this.IXWrap(i,j)
                }
            }
        }

        //console.log(this.ixList);

    }

    createDebugGrid(scene,width,height,blockSize){
        scene.add.grid(width/2,height/2,width, height,blockSize,blockSize,0,255,200);
        
    }


    fluidStep(){
        var N        = this.size;
        var visc     = this.visc;
        var diff     = this.diff;
        var dt       = this.dt;
        var Vx      = this.Vx;
        var Vy      = this.Vy;
        var Vx0     = this.Vx0;
        var Vy0     = this.Vy0;
        var sR       = this.sR;
        var sG       = this.sG;
        var sB       = this.sB;
        var densityR = this.densityR;
        var densityG = this.densityG;
        var densityB = this.densityB;
            
        //this.diffuseVec(1, Vx0, Vx, 2, Vy0, Vy, visc, dt, this.iter, N );
        this.diffuseArr([Vx0,Vy0 ], [Vx, Vy], visc, dt, this.iter, N );
        this.project(Vx0, Vy0, Vx, Vy, this.iter, N);
        
        this.advectArr([Vx,Vy],[Vx0,Vy0], Vx0, Vy0,  dt, N); 
        this.project(Vx, Vy, Vx0, Vy0, this.iter, N);
        
        //process R G B channels of density
        this.diffuseArr([sR,sG,sB], [densityR,densityG,densityB], diff, dt, this.iter, N);
        this.advectArr([densityR,densityG,densityB], [sR,sG,sB], Vx, Vy, dt, N);

  
    }


    diffuseArr(x, x0, diff, dt, iter, N){
        var a = dt * diff * (N - 2) * (N - 2);
        this.lin_solveArr(x, x0, a, 1 + 4 * a, iter, N);
    }

    lin_solve(x, x0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;

        for (var k = 0; k < iter; k++) {     
            for (var j = 0; j < N ; j++) {
                for (var i = 0; i < N; i++) {
                    var idx = this.IXObj(i,j);                    
                    this.tmpar[0][ idx.ix ] = this.calcNeighborVal(idx,x,x0,a, cRecip, N);                    
                }
            }
            for(var i = 0; i < N*N; i++){
                x[i] = this.tmpar[0][i];
            }

        }
    }
    lin_solveArr(x, x0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;
        for (var k = 0; k < iter; k++) {     
            for (var j = 0; j < N ; j++) {
                for (var i = 0; i < N; i++) {
                    var idx = this.IXObj(i,j);
                    for(var ar = 0 ; ar < x.length; ar ++){
                        this.tmpar[ar][ idx.ix ] = this.calcNeighborVal(idx ,x[ar],x0[ar],a, cRecip, N);
                    }
                }
            }      
            for(var i = 0; i < N*N; i++){
                for(var ar = 0 ; ar < x.length; ar ++){
                    x[ar][i] = this.tmpar[ar][i];
                }          
            }
            
        }
    }

    calcNeighborVal(i ,x,x0,a, cRecip, N){
        var ix = i.ix;
        var right = i.right;
        var left = i.left;
        var up = i.up;
        var down = i.down;
        if(!this.solids[ix]){
            return (
                    x0[ix] + a *
                    (   this.solids[right]   ? x[ix] : x[ right ]  + //+ 1 on x axis
                        this.solids[left]    ? x[ix] : x[ left ]  + //-1 on x axis
                        this.solids[down]    ? x[ix] : x[ down ]  + /// +1 on y axis
                        this.solids[up]      ? x[ix] : x[ up ] // -1 on y axis
                    )
                ) * cRecip;
        } else {
            return 0;
        }
            
    }

    project(velocX, velocY, p, div, iter, N)
    {
        var nRecip = 1 / N;   
        for (var j = 0; j < N ; j++) {
            for (var i = 0; i < N; i++) {
                var idx = this.IXObj(i, j);
                div[ idx.ix ] = -0.5*(
                            (  velocX[ idx.right] )
                            -( velocX[ idx.left  ] )
                            +( velocY[ idx.down  ])
                            -( velocY[ idx.up  ] )
                ) * nRecip;
                p[ idx.ix ] = 0;
            }
        }
    
        this.lin_solve( p, div, 1, 6, iter, N);
           
        for (var j = 0; j < N ; j++) {
            for (var i = 0; i < N; i++) {
                var idx = this.IXObj(i, j);
                velocX[ idx.ix ] -= 0.5 * (  p[ idx.right ]- p[ idx.left ]) * N;
                velocY[ idx.ix ] -= 0.5 * (  p[ idx.down ]- p[ idx.up ]) * N;
            }
        }

    }

    advectArr(d, d0, u, v, dt, N)
    {
        var i0, i1, j0, j1;
               
        var s0, s1, t0, t1;
        var  x, y;
        
        var ifloat, jfloat;
        var i, j;
        var dt0 = dt*N;
        var tleft, t, tnext,  vx, vy;

        for(j = 0;  j < N; j++) { 
            for(i = 0; i < N; i++ ) {
                var ix = this.IX(i, j);

                tleft=dt0;
                x=i;y=j;		
        
                while(tleft>0.0000001) {
        
                    //enforce boundry contraints

                    var c = this.checkBounds(x,y);
                    x = c[0]; y = c[1];
                
        
                    i0=Math.floor(x); i1=i0+1;
                    j0=Math.floor(y); j1=j0+1;

                   
                    s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
        
                    vx = -(s0*(t0*u[this.IX(i0,j0)]+t1*u[this.IX(i0,j1)])+
                           s1*(t0*u[this.IX(i1,j0)]+t1*u[this.IX(i1,j1)]));
                            
                    vy = -(s0*(t0*v[this.IX(i0,j0)]+t1*v[this.IX(i0,j1)])+
                           s1*(t0*v[this.IX(i1,j0)]+t1*v[this.IX(i1,j1)]));
        
        
                    var speed2=vx*vx+vy*vy; 
                    if(speed2>0.0000001) tnext=.5/Math.sqrt(speed2);
                    else tnext=tleft;
        
                    t=tnext > tleft ? tleft : tnext;
                    tleft-=t;

                    x+=t*vx;
                    y+=t*vy;
                }
        
        
                var c = this.checkBounds(x,y);
                x = c[0]; y = c[1];
                
        
                i0=Math.floor(x); i1=i0+1;
                j0=Math.floor(y); j1=j0+1;
                s1 = x-i0; s0 = 1-s1; t1 = y-j0; t0 = 1-t1;
                for(var ar = 0 ; ar < d.length; ar ++){
                    d[ar][this.IX(i,j)] =   s0*(t0*d0[ar][this.IX(i0,j0)]+t1*d0[ar][this.IX(i0,j1)])+
                                            s1*(t0*d0[ar][this.IX(i1,j0)]+t1*d0[ar][this.IX(i1,j1)]);            
                }
            }
        }

    }
  

    addItem(x,y,type,amount){
        var ix = this.IXWrap(x,y);
        if(type == "red"){
            this.densityR[ ix ] += amount;
            //this.densityR[ this.IX(i+1,j) ] += amount;
            //this.densityR[ this.IX(i-1,j) ] += amount;
            //this.densityR[ ix + this.size ] += amount;
            //this.densityR[ ix - this.size ] += amount;
        }
        if(type == "green"){
            this.densityG[ ix ] += amount;
            //this.densityG[ this.IX(i+1,j) ] += amount;
            //this.densityG[ this.IX(i-1,j) ] += amount;
            //this.densityG[ ix + this.size ] += amount;
            //this.densityG[ ix - this.size ] += amount;
        }
        if(type == "blue"){
            this.densityB[ ix ] += amount;
            //this.densityB[ this.IX(i+1,j) ] += amount;
            //this.densityB[ this.IX(i-1,j) ] += amount;
            //this.densityB[ ix + this.size ] += amount;
            //this.densityB[ ix - this.size ] += amount;
        }
        if(type == "addwall"){
            this.solids[ ix ] = true;

        }
        if(type == "removewall"){
            this.solids[ ix ] = false;

        }
    }
    addDensity(x, y, amount){
        this.density[ this.IXWrap(x,y) ] += amount;
        
    }
    addVelocity(x, y, amountX, amountY){
        this.Vx[ this.IXWrap(x,y) ] += amountX;
        this.Vy[ this.IXWrap(x,y) ] += amountY;
    }

    getDensity(x,y){       
        return this.densityR[ this.IXWrap(x,y) ]+','+
                    this.densityG[ this.IXWrap(x,y) ]+','+
                    this.densityB[ this.IXWrap(x,y) ];       
    }

    getVelocityX(x,y){
        
        return this.Vx[ this.IXWrap(x,y) ];
        
    }
    getVelocityY(x,y){
        
        return this.Vy[ this.IXWrap(x,y) ];
        
    }

    IX(x, y){
        return this.ixList[x][y].ix;
    }
    IXObj(x, y){        
        return this.ixList[x][y];
    }
    IXWrap(x, y){
        var c = this.checkBounds(x,y);
        return c[0] + c[1] * this.size;;

    }
    checkBounds(x,y){

        if(this.wrap){    
            if(x > this.size-1){
                x = x - this.size; 
            } 
            if(x < 0){
                x = x + this.size;
            }

            if(y > this.size-1){
                y = y - this.size; 
            } 
            if(y < 0){
                y = y + this.size;
            }
            return [x,y];
        }
    }
    render(){
        
        this.graphics.clear();
        

        for(var i = 0; i < this.size; i++){
            for(var j = 0; j < this.size; j++){
                var x = i * this.blockSize;
                var y = j * this.blockSize;

                var ix = this.IX(i,j);
                if( this.solids[ ix ] ){
                    var color = Phaser.Display.Color.GetColor(255,255,255);
                    this.graphics.fillStyle(color);
                    this.graphics.fillRect(x, y, this.blockSize, this.blockSize);
                } else {
                    var r = this.densityR[ ix ];
                    if(r > 255){
                        r = 255;
                    }
                    var g = this.densityG[ ix ];
                    if(g > 255){
                        g = 255;
                    }
                    var b = this.densityB[ ix ];
                    if(b > 255){
                        b = 255;
                    }
                    var color = Phaser.Display.Color.GetColor(r,g,b);
                    this.graphics.fillStyle(color);
                    this.graphics.fillRect(x, y, this.blockSize, this.blockSize);
                }
                
            }    
        }
        
        for(var i = 0; i < this.size; i++){
            for(var j = 0; j < this.size; j++){
                
                var x = i ;
                var y = j ;

                var vx = this.Vx[ this.IX(i,j) ] * 10;
                var vy = this.Vy[ this.IX(i,j) ] * 10;

                var r = Math.abs(vx * 100);
                var b = Math.abs(vy * 100);
                if(r > 255){
                    r = 255;
                }
                
                if(b > 255){
                    b = 255;
                }
                var color = Phaser.Display.Color.GetColor(r,128,b);
                this.graphics.lineStyle(2,color,0.5);
                //if(!(Math.abs(vx) < 0.1 && Math.abs(vy) <= 0.5)){
                    this.graphics.lineBetween(x* this.blockSize, y* this.blockSize, (x+vx)* this.blockSize,(y+vy)* this.blockSize);
                //} 
                
            }    
        }
        
    }

}

