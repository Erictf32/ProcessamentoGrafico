// ========================================
// PROJETO THREE.JS - VISUALIZAÇÃO 3D
// ========================================

// Variáveis globais
let scene, renderer, controls;
let perspectiveCamera, orthographicCamera, currentCamera;
let cube, sphere, plane;
let customShaderMaterial;
let clock;

// Carregar shaders
let vertexShader = '';
let fragmentShader = '';

// ========================================
// INICIALIZAÇÃO
// ========================================

async function init() {
    console.log('Inicializando projeto Three.js...');
    
    // Verificar se Three.js foi carregado
    if (typeof THREE === 'undefined') {
        console.error('Three.js não foi carregado!');
        return;
    }
    
    console.log('Three.js carregado com sucesso, versão:', THREE.REVISION);
    
    // Carregar shaders
    await loadShaders();
    
    // Criar cena
    createScene();
    
    // Criar câmeras
    createCameras();
    
    // Criar renderer
    createRenderer();
    
    // Criar objetos 3D
    createObjects();
    
    // Configurar controles
    setupControls();
    
    // Configurar eventos
    setupEventListeners();
    
    // Iniciar loop de animação
    clock = new THREE.Clock();
    animate();
    
    console.log('Projeto inicializado com sucesso!');
}

// ========================================
// CARREGAMENTO DE SHADERS
// ========================================

async function loadShaders() {
    console.log('Carregando shaders...');
    
    // Primeiro definir shaders fallback
    const fallbackVertexShader = `
        varying vec2 vUv;
        varying vec3 vPosition;
        uniform float time;
        
        void main() {
            vUv = uv;
            vPosition = position;
            vec3 newPosition = position;
            newPosition.z += sin(position.x * 10.0 + time) * 0.1;
            newPosition.z += cos(position.y * 10.0 + time) * 0.1;
            gl_Position = projectionMatrix * modelViewMatrix * vec4(newPosition, 1.0);
        }
    `;
    
    const fallbackFragmentShader = `
        precision mediump float;
        uniform float time;
        varying vec2 vUv;
        varying vec3 vPosition;
        
        void main() {
            float r = sin(vPosition.x * 5.0 + time) * 0.5 + 0.5;
            float g = cos(vPosition.y * 5.0 + time) * 0.5 + 0.5;
            float b = sin(vPosition.z * 5.0 + time * 2.0) * 0.5 + 0.5;
            
            float pattern = sin(vUv.x * 20.0 + time) * cos(vUv.y * 20.0 + time);
            vec3 finalColor = vec3(r, g, b) * (0.7 + pattern * 0.3);
            
            gl_FragColor = vec4(finalColor, 1.0);
        }
    `;

    try {
        console.log('Tentando carregar shaders externos...');
        
        // Carregar vertex shader
        const vertexResponse = await fetch('shaders/custom-vertex.glsl');
        if (vertexResponse.ok) {
            vertexShader = await vertexResponse.text();
            console.log('Vertex shader carregado com sucesso');
        } else {
            throw new Error('Falha ao carregar vertex shader');
        }
        
        // Carregar fragment shader
        const fragmentResponse = await fetch('shaders/custom-fragment.glsl');
        if (fragmentResponse.ok) {
            fragmentShader = await fragmentResponse.text();
            console.log('Fragment shader carregado com sucesso');
        } else {
            throw new Error('Falha ao carregar fragment shader');
        }
        
        console.log('Shaders externos carregados com sucesso');
    } catch (error) {
        console.warn('Erro ao carregar shaders externos:', error.message);
        console.log('Usando shaders fallback...');
        
        // Usar shaders fallback
        vertexShader = fallbackVertexShader;
        fragmentShader = fallbackFragmentShader;
    }
}

// ========================================
// CRIAÇÃO DA CENA
// ========================================

function createScene() {
    console.log('Criando cena...');
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x1a1a1a); // Fundo cinza escuro
    
    // Adicionar iluminação básica
    const ambientLight = new THREE.AmbientLight(0x404040, 0.4);
    scene.add(ambientLight);
    
    const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
    directionalLight.position.set(10, 10, 5);
    scene.add(directionalLight);
    
    console.log('Cena criada com sucesso');
}

// ========================================
// CRIAÇÃO DAS CÂMERAS
// ========================================

function createCameras() {
    console.log('Criando câmeras...');
    const aspect = window.innerWidth / window.innerHeight;
    
    // Câmera perspectiva
    perspectiveCamera = new THREE.PerspectiveCamera(75, aspect, 0.1, 1000);
    perspectiveCamera.position.set(5, 5, 5);
    perspectiveCamera.lookAt(0, 0, 0);
    
    // Câmera ortográfica
    const frustumSize = 10;
    orthographicCamera = new THREE.OrthographicCamera(
        frustumSize * aspect / -2,
        frustumSize * aspect / 2,
        frustumSize / 2,
        frustumSize / -2,
        1,
        1000
    );
    orthographicCamera.position.set(5, 5, 5);
    orthographicCamera.lookAt(0, 0, 0);
    
    // Câmera atual (inicia com perspectiva)
    currentCamera = perspectiveCamera;
    
    console.log('Câmeras criadas com sucesso');
}

// ========================================
// CRIAÇÃO DO RENDERER
// ========================================

function createRenderer() {
    console.log('Criando renderer...');
    renderer = new THREE.WebGLRenderer({ 
        antialias: true,
        alpha: true 
    });
    renderer.setSize(window.innerWidth, window.innerHeight);
    renderer.setPixelRatio(window.devicePixelRatio);
    
    // Adicionar canvas ao DOM
    const container = document.getElementById('canvas-container');
    container.appendChild(renderer.domElement);
    
    console.log('Renderer criado com sucesso');
}

// ========================================
// CRIAÇÃO DOS OBJETOS 3D
// ========================================

function createObjects() {
    console.log('Criando objetos 3D...');
    
    // 1. CUBO COM SHADER CUSTOMIZADO
    createCubeWithCustomShader();
    
    // 2. ESFERA COM TEXTURA
    createSphereWithTexture();
    
    // 3. PLANO COM TEXTURA XADREZ
    createPlaneWithCheckerboard();
    
    console.log('Todos os objetos 3D criados com sucesso');
}

function createCubeWithCustomShader() {
    console.log('Criando cubo com shader customizado...');
    
    try {
        // Geometria do cubo
        const geometry = new THREE.BoxGeometry(2, 2, 2);
        console.log('Geometria do cubo criada');
        
        // Material com shader customizado (usando ShaderMaterial em vez de RawShaderMaterial)
        customShaderMaterial = new THREE.ShaderMaterial({
            vertexShader: vertexShader,
            fragmentShader: fragmentShader,
            uniforms: {
                time: { value: 0.0 }
            }
        });
        console.log('Material shader criado');
        
        // Criar cubo
        cube = new THREE.Mesh(geometry, customShaderMaterial);
        cube.position.set(-3, 0, 0);
        cube.scale.set(1.2, 1.2, 1.2);
        scene.add(cube);
        
        console.log('Cubo criado com sucesso! Posição:', cube.position);
    } catch (error) {
        console.error('Erro ao criar cubo com shader:', error);
        
        // Fallback: criar cubo com material básico
        console.log('Criando cubo com material básico como fallback...');
        const geometry = new THREE.BoxGeometry(2, 2, 2);
        const material = new THREE.MeshLambertMaterial({ color: 0xff6b6b });
        cube = new THREE.Mesh(geometry, material);
        cube.position.set(-3, 0, 0);
        cube.scale.set(1.2, 1.2, 1.2);
        scene.add(cube);
        console.log('Cubo fallback criado com sucesso!');
    }
}

function createSphereWithTexture() {
    console.log('Criando esfera...');
    
    try {
        // Geometria da esfera
        const geometry = new THREE.SphereGeometry(1.5, 32, 32);
        
        // Criar textura gradiente
        if (typeof createGradientTexture === 'function') {
            const canvas = createGradientTexture(256, 256);
            const texture = new THREE.CanvasTexture(canvas);
            const material = new THREE.MeshLambertMaterial({ map: texture });
            
            // Criar esfera
            sphere = new THREE.Mesh(geometry, material);
        } else {
            // Fallback sem textura
            const material = new THREE.MeshLambertMaterial({ color: 0x4ecdc4 });
            sphere = new THREE.Mesh(geometry, material);
        }
        
        sphere.position.set(3, 0, 0);
        sphere.scale.set(1.0, 1.0, 1.0);
        scene.add(sphere);
        
        console.log('Esfera criada com sucesso');
    } catch (error) {
        console.error('Erro ao criar esfera:', error);
    }
}

function createPlaneWithCheckerboard() {
    console.log('Criando plano...');
    
    try {
        // Geometria do plano
        const geometry = new THREE.PlaneGeometry(32, 32);
        
        // Criar textura xadrez
        if (typeof createCheckerboardTexture === 'function') {
            const canvas = createCheckerboardTexture(512, 512);
            const texture = new THREE.CanvasTexture(canvas);
            texture.wrapS = THREE.RepeatWrapping;
            texture.wrapT = THREE.RepeatWrapping;
            texture.repeat.set(2, 2);
            const material = new THREE.MeshLambertMaterial({ map: texture });
            
            // Criar plano
            plane = new THREE.Mesh(geometry, material);
        } else {
            // Fallback sem textura
            const material = new THREE.MeshLambertMaterial({ color: 0x888888 });
            plane = new THREE.Mesh(geometry, material);
        }
        
        plane.rotation.x = -Math.PI / 2;
        plane.position.set(0, -2, 0);
        scene.add(plane);
        
        console.log('Plano criado com sucesso');
    } catch (error) {
        console.error('Erro ao criar plano:', error);
    }
}

// ========================================
// CONFIGURAÇÃO DOS CONTROLES
// ========================================

function setupControls() {
    console.log('Configurando controles...');
    
    try {
        // Verificar se OrbitControls está disponível
        if (typeof THREE.OrbitControls !== 'undefined') {
            controls = new THREE.OrbitControls(currentCamera, renderer.domElement);
        } else if (typeof OrbitControls !== 'undefined') {
            controls = new OrbitControls(currentCamera, renderer.domElement);
        } else {
            console.warn('OrbitControls não disponível');
            return;
        }
        
        controls.enableDamping = true;
        controls.dampingFactor = 0.05;
        controls.screenSpacePanning = false;
        controls.minDistance = 3;
        controls.maxDistance = 50;
        controls.maxPolarAngle = Math.PI / 2;
        
        console.log('Controles configurados com sucesso');
    } catch (error) {
        console.error('Erro ao configurar controles:', error);
    }
}

// ========================================
// EVENTOS
// ========================================

function setupEventListeners() {
    console.log('Configurando event listeners...');
    
    // Redimensionamento da janela
    window.addEventListener('resize', onWindowResize);
    
    // Troca de câmera com tecla C
    window.addEventListener('keydown', onKeyDown);
    
    console.log('Event listeners configurados');
}

function onWindowResize() {
    const aspect = window.innerWidth / window.innerHeight;
    
    // Atualizar câmera perspectiva
    perspectiveCamera.aspect = aspect;
    perspectiveCamera.updateProjectionMatrix();
    
    // Atualizar câmera ortográfica
    const frustumSize = 10;
    orthographicCamera.left = frustumSize * aspect / -2;
    orthographicCamera.right = frustumSize * aspect / 2;
    orthographicCamera.top = frustumSize / 2;
    orthographicCamera.bottom = frustumSize / -2;
    orthographicCamera.updateProjectionMatrix();
    
    // Atualizar renderer
    renderer.setSize(window.innerWidth, window.innerHeight);
}

function onKeyDown(event) {
    switch (event.code) {
        case 'KeyC':
            toggleCamera();
            break;
    }
}

function toggleCamera() {
    if (currentCamera === perspectiveCamera) {
        currentCamera = orthographicCamera;
        console.log('Câmera alterada para: Ortográfica');
    } else {
        currentCamera = perspectiveCamera;
        console.log('Câmera alterada para: Perspectiva');
    }
    
    // Atualizar controles
    if (controls) {
        controls.object = currentCamera;
        controls.update();
    }
}

// ========================================
// LOOP DE ANIMAÇÃO
// ========================================

function animate() {
    requestAnimationFrame(animate);
    
    const elapsedTime = clock.getElapsedTime();
    
    // Animações dos objetos
    updateAnimations(elapsedTime);
    
    // Atualizar controles
    if (controls) {
        controls.update();
    }
    
    // Renderizar cena
    renderer.render(scene, currentCamera);
}

function updateAnimations(time) {
    // Animação do cubo (rotação e shader)
    if (cube) {
        cube.rotation.x = time * 0.5;
        cube.rotation.y = time * 0.7;
        
        // Atualizar uniform do shader se disponível
        if (customShaderMaterial && customShaderMaterial.uniforms && customShaderMaterial.uniforms.time) {
            customShaderMaterial.uniforms.time.value = time;
        }
    }
    
    // Animação da esfera (movimento de quique e rotação)
    if (sphere) {
        // Movimento de quique: a esfera "bate" no plano sem atravessá-lo
        // Plano está em y = -2, esfera tem raio 1.5
        // Centro da esfera deve ficar entre y = -0.5 (tocando o plano) e y = 1.0 (altura do quique)
        const bounceHeight = 1.5; // Altura máxima do quique
        const minHeight = -0.5;   // Altura mínima (tocando o plano)
        
        // Usar abs(sin) para criar movimento de quique mais realista
        const bounceMotion = Math.abs(Math.sin(time * 3));
        sphere.position.y = minHeight + (bounceMotion * bounceHeight);
        
        sphere.rotation.y = time;
    }
    
    // Animação do plano (rotação suave)
    if (plane) {
        plane.rotation.z = Math.sin(time * 0.3) * 0.1;
    }
}

// ========================================
// INICIALIZAÇÃO DO PROJETO
// ========================================

// Iniciar quando a página carregar
window.addEventListener('DOMContentLoaded', init); 