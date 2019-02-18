//https://mikeash.com/pyblog/fluid-simulation-for-dummies.html
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
        
        this.iter = 16;
        //size of the grid
        this.size = gridWidth/blockSize;
        this.gridsize = this.size*this.size;
        this.N = this.size;
        //the timestamp, diffusion and viscosity of this fluidgrid
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        this.s = new Array(this.gridsize).fill(0);
        this.density = new Array(this.gridsize).fill(0);
        this.Vx = new Array(this.gridsize).fill(0);
        this.Vx0 = new Array(this.gridsize).fill(0);
        this.Vy = new Array(this.gridsize).fill(0);
        this.Vy0 = new Array(this.gridsize).fill(0);

        this.graphics = this.scene.add.graphics();
        //grid size is size * size
        //this.grid = this.initaliseFluidGrid(this.size);
        //this.densityVisual = new Array(this.gridsize).fill( new Phaser.Geom.Rectangle(0,0,this.blockSize,this.blockSize)); 
        //this.vectorVisual = new Array(this.gridsize).fill( new Phaser.Geom.Line(0,0,0,0)); 
    }

    createDebugGrid(scene,width,height,blockSize){
        scene.add.grid(width/2,height/2,width, height,blockSize,blockSize,0,255,200);
        
    }

    initaliseFluidGrid(size){
        //initialise our grid in a 1D array
        var grid = [];

        for(var i = 0; i < size * size; i++){
            var Vx = new Vector(0,0);
            grid[i] = {
                s: null,
                density: null,
                Vx: null,
                Vx0: null,
                Vy: null,
                Vy0: null
                
            }
        }
        return grid;
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
        var s       = this.s;
        var density = this.density;
    
        this.diffuse(1, Vx0, Vx, visc, dt, this.iter, N); // x
        this.diffuse(2, Vy0, Vy, visc, dt, this.iter, N); // y
        
        
        this.project(Vx0, Vy0, Vx, Vy, this.iter, N);
        
        this.advect(1, Vx, Vx0, Vx0, Vy0,  dt, N); // x
        this.advect(2, Vy, Vy0, Vx0, Vy0,  dt, N); // y
        
        
        this.project(Vx, Vy, Vx0, Vy0, this.iter, N);
        
        this.diffuse(0, s, density, diff, dt, this.iter, N);
        this.advect(0, density, s, Vx, Vy, dt, N);
    }

    //type, oldvel, newvel, visc, timestep, iteration, size
    diffuse(b, x, x0, diff, dt, iter, N){
        var a = dt * diff * (N - 2) * (N - 2);
        this.lin_solve(b, x, x0, a, 1 + 6 * a, iter, N);
    }



    lin_solve(b, x, x0,  a,  c,  iter,  N){
        var cRecip = 1.0 / c;
        for (var k = 0; k < iter; k++) {     
            for (var j = 1; j < N - 1; j++) {
                for (var i = 1; i < N - 1; i++) {
                    x[this.IX(i, j)] =
                        (x0[this.IX(i, j)]
                            + a*(    x[this.IX(i+1, j   )]
                                    +x[this.IX(i-1, j   )]
                                    +x[this.IX(i  , j+1 )]
                                    +x[this.IX(i  , j-1 )]
                                    //+x[this.IX(i  , j   )] //do i need? would have been z -1
                        )) * cRecip;
                }
            }            
            this.set_bnd(b, x, N);
        }
    }

    project(velocX, velocY, p, div, iter, N)
    {
        
        for (var j = 1; j < N - 1; j++) {
            for (var i = 1; i < N - 1; i++) {
                div[this.IX(i, j )] = -0.5*(
                            velocX[this.IX(i+1, j    )]
                            -velocX[this.IX(i-1, j    )]
                            +velocY[this.IX(i  , j+1  )]
                            -velocY[this.IX(i  , j-1  )]
                    )/N;
                p[this.IX(i, j )] = 0;
            }
        }
        
        this.set_bnd(0, div, N); 
        this.set_bnd(0, p, N);
        this.lin_solve(0, p, div, 1, 6, iter, N);
        
        
        for (var j = 1; j < N - 1; j++) {
            for (var i = 1; i < N - 1; i++) {
                velocX[this.IX(i, j)] -= 0.5 * (  p[this.IX(i+1, j )]
                                                -p[this.IX(i-1, j)]) * N;
                velocY[this.IX(i, j)] -= 0.5 * (  p[this.IX(i, j+1)]
                                                -p[this.IX(i, j-1)]) * N;
            }
        }
        
        this.set_bnd(1, velocX, N);
        this.set_bnd(2, velocY, N);
    }

    advect(b, d, d0,  velocX, velocY, dt, N)
    {
        var i0, i1, j0, j1;
        
        var dtx = dt * (N - 2);
        var dty = dt * (N - 2);
        
        var s0, s1, t0, t1;
        var tmp1, tmp2, x, y;
        
        var Nfloat = N;
        var ifloat, jfloat;
        var i, j;
        

        for(j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
            for(i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[this.IX(i, j)];
                tmp2 = dty * velocY[this.IX(i, j)];

                x    = ifloat - tmp1; 
                y    = jfloat - tmp2;
                
                if(x < 0.5) x = 0.5; 
                if(x > Nfloat + 0.5) x = Nfloat + 0.5; 
                i0 = Math.floor(x); 
                i1 = i0 + 1.0;
                if(y < 0.5) y = 0.5; 
                if(y > Nfloat + 0.5) y = Nfloat + 0.5; 
                j0 = Math.floor(y);
                j1 = j0 + 1.0; 
                
                s1 = x - i0; 
                s0 = 1.0 - s1; 
                t1 = y - j0; 
                t0 = 1.0 - t1;
                
                var i0i = Math.floor(i0);
                var i1i = Math.floor(i1);
                var j0i = Math.floor(j0);
                var j1i = Math.floor(j1);
                
                d[this.IX(i, j)] = 
                
                    s0 * ( t0 * d0[this.IX(i0i, j0i)]  +  t1 * d0[this.IX(i0i, j1i)] ) +
                    s1 * ( t0 * d0[this.IX(i1i, j0i)]  +  t1 * d0[this.IX(i1i, j1i)] );
            }
        }
    
        this.set_bnd(b, d, N);
    }

    set_bnd(b, x, N){

        
        for(var i = 1; i < N - 1; i++) {
            x[this.IX(i, 0  )] = b == 2 ? -x[this.IX(i, 1  )] : x[this.IX(i, 1  )];
            x[this.IX(i, N-1)] = b == 2 ? -x[this.IX(i, N-2)] : x[this.IX(i, N-2)];
        }
    
    
        for(var j = 1; j < N - 1; j++) {
            x[this.IX(0  , j)] = b == 1 ? -x[this.IX(1  , j)] : x[this.IX(1  , j)];
            x[this.IX(N-1, j)] = b == 1 ? -x[this.IX(N-2, j)] : x[this.IX(N-2, j)];
        }
        
        
        x[this.IX(0, 0)]     = 0.5 * (x[this.IX(1, 0)]     + x[this.IX(0, 1)]     + x[this.IX(0, 0)]);
        x[this.IX(0, N-1)]   = 0.5 * (x[this.IX(1, N-1)]   + x[this.IX(0, N-2)]   + x[this.IX(0, N-1)]);
        x[this.IX(N-1, 0)]   = 0.5 * (x[this.IX(N-2, 0)]   + x[this.IX(N-1, 1)]   + x[this.IX(N-1, 0)]);
        x[this.IX(N-1, N-1)] = 0.5 * (x[this.IX(N-2, N-1)] + x[this.IX(N-1, N-2)] + x[this.IX(N-1, N-1)]);
    }

    addDensity(x, y, amount){
        this.density[ this.IX(x,y) ] += amount;
        
    }
    addVelocity(x, y, amountX, amountY){
        this.Vx[ this.IX(x,y) ] += amountX;
        this.Vy[ this.IX(x,y) ] += amountY;
    }

    getDensity(x,y){
        return this.density[ this.IX(x,y) ];
    }

    getVelocityX(x,y){
        return this.Vx[ this.IX(x,y) ];
    }
    getVelocityY(x,y){
        return this.Vy[ this.IX(x,y) ];
    }

    IX(x, y){
        //converts the X, Y coords to a 1D array, size being the length of Y

        return Math.max(0, Math.min( x, this.size-1)) + Math.max(0, Math.min( y, this.size-1)) * this.size;
    }

    
    render(){
        
        this.graphics.clear();
        

        for(var i = 0; i < this.size; i++){
            for(var j = 0; j < this.size; j++){
                
                var x = i * this.blockSize;
                var y = j * this.blockSize;
                
                var d = this.density[ this.IX(i,j) ];
                if(d > 255){
                    d = 255;
                }
                var color = Phaser.Display.Color.GetColor(d,0,0);
                this.graphics.fillStyle(color);
                this.graphics.fillRect(x, y, this.blockSize, this.blockSize);
                
            }    
        }
        
        for(var i = 0; i < this.size; i++){
            for(var j = 0; j < this.size; j++){
                
                var x = i ;
                var y = j ;

                var vx = this.Vx[ this.IX(i,j) ];
                var vy = this.Vy[ this.IX(i,j) ];

                this.graphics.lineStyle(2,0x44ff22,0.1);
                if(!(Math.abs(vx) < 0.1 && Math.abs(vy) <= 0.1)){
                    this.graphics.lineBetween(x* this.blockSize, y* this.blockSize, (x+vx)* this.blockSize,(y+vy)* this.blockSize);
                } 
                
            }    
        }
        
    }

}

