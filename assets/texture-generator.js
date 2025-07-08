// Gerador de Textura Procedural
function createCheckerboardTexture(width = 256, height = 256) {
    const canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;
    
    const context = canvas.getContext('2d');
    
    // Criar padrão xadrez
    const tileSize = 32;
    for (let x = 0; x < width; x += tileSize) {
        for (let y = 0; y < height; y += tileSize) {
            const isEven = ((x / tileSize) + (y / tileSize)) % 2 === 0;
            context.fillStyle = isEven ? '#ff6b6b' : '#4ecdc4';
            context.fillRect(x, y, tileSize, tileSize);
        }
    }
    
    return canvas;
}

function createGradientTexture(width = 256, height = 256) {
    const canvas = document.createElement('canvas');
    canvas.width = width;
    canvas.height = height;
    
    const context = canvas.getContext('2d');
    
    // Criar gradiente radial
    const gradient = context.createRadialGradient(
        width/2, height/2, 0,
        width/2, height/2, width/2
    );
    gradient.addColorStop(0, '#ffeb3b');
    gradient.addColorStop(0.5, '#ff9800');
    gradient.addColorStop(1, '#e91e63');
    
    context.fillStyle = gradient;
    context.fillRect(0, 0, width, height);
    
    return canvas;
} 